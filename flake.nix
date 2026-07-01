{
  description = "Optional reproducible development shell for cancer_research";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    rust-overlay.url = "github:oxalica/rust-overlay";
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
    rust-overlay,
  }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        overlays = [(import rust-overlay)];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
        python = pkgs.python312;
        rustToolchain = pkgs.rust-bin.fromRustupToolchainFile ./simulations/rust-toolchain.toml;
        tex = pkgs.texlive.combine {
          inherit (pkgs.texlive)
            scheme-medium
            latexmk
            collection-fontsrecommended
            collection-latexextra;
        };
      in {
        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            python
            uv
            rustToolchain
            maturin
            pkg-config
            openssl
            zlib
            libxml2
            libxslt
            graphviz
            git
            git-lfs
            tex
          ];

          shellHook = ''
            export UV_PYTHON_DOWNLOADS=never
            export UV_PYTHON="${python}/bin/python3"
            export PYO3_PYTHON="${python}/bin/python3"

            echo "Entered cancer_research nix dev shell"
            echo "Run: uv sync --frozen"
          '';
        };
      });
}
