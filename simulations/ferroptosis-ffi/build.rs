use std::path::PathBuf;

fn main() {
    let crate_dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let config_path = crate_dir.join("cbindgen.toml");
    let output_path = crate_dir.join("include").join("ferroptosis.h");

    // Ensure include directory exists
    std::fs::create_dir_all(crate_dir.join("include")).ok();

    let config = cbindgen::Config::from_file(&config_path)
        .unwrap_or_default();

    match cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(config)
        .generate()
    {
        Ok(bindings) => {
            bindings.write_to_file(&output_path);
        }
        Err(e) => {
            eprintln!("cargo:warning=cbindgen failed: {e}");
        }
    }
}
