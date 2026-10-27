fn main() -> std::process::ExitCode {
    let cli = rsomics_help::parse::<rsomics_rnaseq_qc::cli::Cli>();
    cli.run()
}
