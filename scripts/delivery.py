"""GitHub source/native delivery; never a Cargo registry package."""

import argparse
import base64
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tarfile

PRODUCT = 'rsomics-rnaseq-qc'
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METRICS = ['mapping', 'strandedness', 'distribution', 'coverage', 'integrity']


def identity(path):
    body = Path(path).read_bytes()
    return dict(bytes=len(body), sha256=hashlib.sha256(body).hexdigest())


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def canonical(value):
    return (json.dumps(value, indent=2, sort_keys=True) + '\n').encode()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def source_manifest():
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    require(commit == os.environ['GITHUB_SHA'], 'source commit differs from the actual CI head')
    metadata = json.loads(subprocess.check_output(['cargo', 'metadata', '--locked', '--no-deps', '--format-version', '1'], cwd=ROOT))
    package, = [p for p in metadata['packages'] if p['name'] == PRODUCT]
    records = subprocess.check_output(['git', 'ls-tree', '-rz', '--full-tree', commit], cwd=ROOT)
    files = {}
    for record in records.rstrip(b'\0').split(b'\0'):
        header, name = record.split(b'\t', 1)
        mode, kind, blob = header.decode().split()
        name = name.decode()
        require(kind == 'blob' and mode in ('100644', '100755'), 'unsupported Git source member: ' + name)
        body = (ROOT / name).read_bytes()
        observed_blob = hashlib.sha1(b'blob ' + str(len(body)).encode() + b'\0' + body).hexdigest()
        require(observed_blob == blob, 'worktree body differs from exact Git source: ' + name)
        require(bool((ROOT / name).stat().st_mode & 0o111) == (mode == '100755'), 'source mode differs: ' + name)
        files[name] = dict(bytes=len(body), sha256=hashlib.sha256(body).hexdigest(), git_blob=blob, git_mode=mode)
    return dict(schema='rsomics.github-source.v1', product=PRODUCT, version=package['version'],
        prefix=PRODUCT + '-' + package['version'] + '-source', source_commit=commit,
        source_tree=subprocess.check_output(['git', 'rev-parse', commit + '^{tree}'], cwd=ROOT, text=True).strip(),
        files=files, registry_package=False, dependencies_vendored=False,
        dependency_policy='Original Cargo.toml/Cargo.lock preserve public pinned Git dependencies; fresh builds need network/cache access.')


def write_source_tar(root, manifest, output):
    root, output = Path(root), Path(output)
    bodies = {}
    for name, row in manifest['files'].items():
        body = (root / name).read_bytes()
        require(dict(bytes=len(body), sha256=hashlib.sha256(body).hexdigest()) == {k: row[k] for k in ('bytes', 'sha256')},
                'original source body changed: ' + name)
        bodies[name] = body
    bodies['SOURCE.json'] = canonical(manifest)
    directories = {manifest['prefix']}
    for name in bodies:
        path = PurePosixPath(manifest['prefix']) / name
        require(not path.is_absolute() and '..' not in path.parts, 'unsafe source path')
        directories.update(str(p) for p in path.parents if str(p) != '.')
    with output.open('xb') as stream, tarfile.open(fileobj=stream, mode='w', format=tarfile.PAX_FORMAT) as archive:
        for name in sorted(directories):
            row = tarfile.TarInfo(name)
            row.type, row.mode, row.mtime = tarfile.DIRTYPE, 0o755, 0
            archive.addfile(row)
        for name, body in sorted(bodies.items()):
            row = tarfile.TarInfo(manifest['prefix'] + '/' + name)
            row.size, row.mtime = len(body), 0
            row.mode = int(manifest['files'][name]['git_mode'], 8) & 0o777 if name != 'SOURCE.json' else 0o644
            archive.addfile(row, io.BytesIO(body))


def extract_source_tar(archive_path, manifest, destination):
    destination = Path(destination)
    destination.mkdir()
    with tarfile.open(archive_path, 'r:') as archive:
        members = archive.getmembers()
        require(len(members) == len({m.name for m in members}), 'duplicate source TAR member')
        original = manifest['files']
        expected = {manifest['prefix'] + '/' + name for name in [*original, 'SOURCE.json']}
        require({m.name for m in members if m.isfile()} == expected, 'source TAR original membership differs')
        require(archive.extractfile(manifest['prefix'] + '/SOURCE.json').read() == canonical(manifest), 'source provenance differs')
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and '..' not in path.parts and path.parts[0] == manifest['prefix'], 'unsafe source TAR path')
            require(member.isfile() or member.isdir(), 'unsupported source TAR type')
            require(member.uid == member.gid == member.mtime == 0 and member.uname == member.gname == '', 'source headers are not deterministic')
            target = destination / member.name
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                target.chmod(0o755)
            else:
                body = archive.extractfile(member).read()
                name = '/'.join(path.parts[1:])
                if name in original:
                    require(dict(bytes=len(body), sha256=hashlib.sha256(body).hexdigest()) == {k: original[name][k] for k in ('bytes', 'sha256')},
                            'extracted original body differs: ' + name)
                    require(member.mode == int(original[name]['git_mode'], 8) & 0o777, 'extracted original mode differs: ' + name)
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open('xb') as stream:
                    stream.write(body)
                target.chmod(member.mode)
    return destination / manifest['prefix']


def source():
    temporary = Path(os.environ['RUNNER_TEMP'])
    delivery = temporary / 'delivery'
    delivery.mkdir()
    metadata_body = subprocess.check_output(['cargo', 'metadata', '--locked', '--format-version', '1'], cwd=ROOT)
    tree_body = subprocess.check_output(['cargo', 'tree', '--locked', '--target', os.environ['EXPECTED_HOST'],
        '--edges', 'normal,build', '--prefix', 'none', '--format', '{p}'], cwd=ROOT)
    (delivery / 'native-dependencies.metadata.json').write_bytes(metadata_body)
    (delivery / 'native-dependencies.tree.txt').write_bytes(tree_body)
    metadata = json.loads(metadata_body)
    by_name = {(p['name'], p['version']): p for p in metadata['packages']}
    enabled = set()
    for line in tree_body.decode().splitlines():
        name, version = re.match(r'([^ ]+) v([^ ]+)', line).groups()
        if name != PRODUCT:
            enabled.add(by_name[(name, version)]['id'])
    inventory = json.loads((ROOT / 'LICENSES/provenance.json').read_bytes())
    expected = {key for key, row in inventory['packages'].items() if os.environ['EXPECTED_HOST'] in row['targets']}
    require(enabled == expected,
        'actual native enabled dependency graph differs from the notice inventory')
    require(identity(ROOT / 'Cargo.lock')['sha256'] == inventory['lock_sha256'], 'notice inventory lock differs')
    save(delivery / 'native-dependencies.json', dict(target=os.environ['EXPECTED_HOST'], enabled_packages=sorted(enabled),
        metadata_sha256=hashlib.sha256(metadata_body).hexdigest(), tree_sha256=hashlib.sha256(tree_body).hexdigest(),
        inventory_sha256=identity(ROOT / 'LICENSES/provenance.json')['sha256'],
        actual_native_enabled_graph_verified=True, linked_binary_inventory=False))
    manifest = source_manifest()
    archive = delivery / (manifest['prefix'] + '.tar')
    write_source_tar(ROOT, manifest, archive)
    extracted = extract_source_tar(archive, manifest, temporary / 'source')
    save(temporary / 'source-verification.json', dict(**identity(archive), name=archive.name,
        original_files=len(manifest['files']), source_commit=manifest['source_commit'], source_tree=manifest['source_tree'],
        all_original_bodies_modes_verified=True, deterministic_headers=True,
        registry_package=False, dependencies_vendored=False))
    with Path(os.environ['GITHUB_ENV']).open('a') as stream:
        stream.write('SOURCE_DIRECTORY=' + str(extracted) + '\n')
        stream.write('SOURCE_TAR=' + str(archive) + '\n')
        stream.write('PACKAGE_VERSION=' + manifest['version'] + '\n')


def decode_fixture(path, output):
    model = json.loads(Path(path).read_bytes())
    require(model['schema'] == 'rsomics.rnaseq-qc.delivery-fixture.v1', 'delivery fixture schema differs')
    output = Path(output)
    output.mkdir()
    for name, row in model['files'].items():
        require(name in ('paired.bam', 'genes.bed12'), 'unexpected delivery fixture member')
        body = base64.b64decode(row['base64'], validate=True)
        require(dict(bytes=len(body), sha256=hashlib.sha256(body).hexdigest()) == {k: row[k] for k in ('bytes', 'sha256')},
                'delivery fixture body differs: ' + name)
        with (output / name).open('xb') as stream:
            stream.write(body)
    save(output / 'PROVENANCE.json', {k:v for k,v in model.items() if k != 'files'})
    return model


def check_output(directory, metrics, records, workers):
    directory = Path(directory)
    manifest = json.loads((directory / 'manifest.json').read_bytes())
    report = json.loads((directory / 'report.json').read_bytes())
    require(manifest['status'] == 'complete' and manifest['bam_traversals'] == 1
        and manifest['original_records'] == records and manifest['inflate_workers'] == workers,
        'installed product input/worker/transaction contract differs')
    require(manifest['selected_metrics'] == report['metrics'] == metrics, 'installed metric selection differs')
    require(all(name in report for name in metrics), 'installed report is missing a selected metric')
    expected = {row['path'] for row in manifest['outputs']} | {'manifest.json'}
    require({p.name for p in directory.iterdir()} == expected, 'installed output set differs from its manifest')
    for row in manifest['outputs']:
        require(identity(directory / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}, 'installed output checksum differs')
    return report, manifest


def smoke(binary):
    binary = Path(binary).resolve(strict=True)
    temporary = Path(os.environ['RUNNER_TEMP'])
    directory = temporary / 'delivery/lifecycle'
    directory.mkdir(parents=True)
    shutil.copy2(binary, directory / 'installed-binary')
    save(directory / 'installed-binary.json', identity(binary))
    fixture = decode_fixture(ROOT / 'tests/fixtures/delivery/paired.json', directory / 'inputs')
    junction = ROOT / 'tests/fixtures/junctions/informative-paired'
    calls = [('--version',), ('--help',),
        ('report', '--bam', str(directory / 'inputs/paired.bam'), '--genes', str(directory / 'inputs/genes.bed12'),
         '--out', str(directory / 'default-five'), '--threads', '4'),
        ('junctions', '--bam', str(junction / 'input.bam'), '--genes', str(junction / 'genes.bed12'),
         '--out', str(directory / 'junctions'), '--profile', 'rseqc-5.0.5', '--threads', '4'),
        ('report', '--bam', str(junction / 'input.bam'), '--genes', str(junction / 'genes.bed12'),
         '--out', str(directory / 'selected'), '--profile', 'rseqc-5.0.5', '--threads', '4', '--metrics', 'mapping,junctions')]
    records = []
    for index, argv in enumerate(calls):
        command = [str(binary), *argv]
        with (directory / f'{index}.stdout').open('xb') as out, (directory / f'{index}.stderr').open('xb') as err:
            result = subprocess.run(command, stdout=out, stderr=err, cwd=directory)
        records.append(dict(argv=command, returncode=result.returncode))
        save(directory / 'commands.json', records)
        require(result.returncode == 0, 'installed command failed: ' + str(command))
    default, _ = check_output(directory / 'default-five', DEFAULT_METRICS, fixture['original_records'], 4)
    require('junctions' not in default, 'default report unexpectedly includes junctions')
    focused, _ = check_output(directory / 'junctions', ['junctions'], 30, 4)
    selected, _ = check_output(directory / 'selected', ['mapping', 'junctions'], 30, 4)
    expected = json.loads((junction / 'observation.json').read_bytes())
    require(focused['junctions']['rows'] == expected['rows'] and selected['junctions'] == focused['junctions'],
        'installed ordered junction rows differ from captured fixture/shared-scan selection')
    require([focused['junctions']['counts'][k] for k in ('total_events', 'known_events', 'partial_novel_events', 'novel_events', 'filtered_events')]
        == [22, 12, 2, 6, 2], 'installed junction counters differ')
    save(directory / 'VERIFIED.json', dict(installed_binary=identity(binary), calls=records,
        default_metrics=DEFAULT_METRICS, default_records=175, junction_records=30,
        original_outputs_retained=True, scope='installed source-TAR lifecycle; not a new source or performance campaign'))


def bundle(binary):
    temporary = Path(os.environ['RUNNER_TEMP'])
    binary = Path(binary).resolve(strict=True)
    version, target = os.environ['PACKAGE_VERSION'], os.environ['NATIVE_TARGET']
    delivery = temporary / 'delivery'
    stem = PRODUCT + '-' + version + '-' + target
    folder = delivery / stem
    folder.mkdir()
    shutil.copy2(binary, folder / PRODUCT)
    for name in ('README.md', 'MERGE_RECORD.md', 'LICENSE-MIT', 'LICENSE-APACHE-2.0', 'THIRD_PARTY_NOTICES.md'):
        shutil.copy2(ROOT / name, folder / name)
    shutil.copytree(ROOT / 'LICENSES', folder / 'LICENSES')
    provenance = dict(product=PRODUCT, version=version, target=target, head_sha=os.environ['GITHUB_SHA'],
        run_id=int(os.environ['GITHUB_RUN_ID']), run_attempt=int(os.environ['GITHUB_RUN_ATTEMPT']),
        native_host=os.environ['EXPECTED_HOST'], rustc=subprocess.check_output(['rustc', '-vV'], text=True),
        source=json.loads((temporary / 'source-verification.json').read_bytes()),
        binary=identity(binary), lock_sha256=identity(ROOT / 'Cargo.lock')['sha256'],
        native_dependency_graph=json.loads((delivery / 'native-dependencies.json').read_bytes()),
        historical_default_five_performance_run=40592936910,
        historical_default_five_binary_sha256='48091abdb22d5d91e18a733848798b8a9ae3dd83d3d8117dd77feca08fd13707',
        historical_junction_performance_run=41291319388,
        historical_junction_binary_sha256='eb472af80dfa26a00819448f895ccbcceb570de3ad629ac81f9507c7534282ae',
        this_binary_remeasured=False, registry_package=False, signed_release=False)
    save(folder / 'BUILD.json', provenance)
    for path in [folder, *folder.rglob('*')]:
        path.chmod(0o755 if path.is_dir() or path.name == PRODUCT else 0o644)
    archive = delivery / (stem + '.tar.gz')
    with tarfile.open(archive, 'w:gz') as output:
        output.add(folder, arcname=stem)
    source_tar = Path(os.environ['SOURCE_TAR'])
    (delivery / 'SHA256SUMS').write_text(''.join(identity(p)['sha256'] + '  ' + p.name + '\n' for p in (archive, source_tar)))
    print(json.dumps(provenance, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('source', 'fixture', 'smoke', 'bundle'))
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if args.action == 'source':
        source()
    elif args.action == 'fixture':
        decode_fixture(ROOT / 'tests/fixtures/delivery/paired.json', args.out)
    elif args.action == 'smoke':
        smoke(args.binary)
    else:
        bundle(args.binary)


if __name__ == '__main__':
    main()
