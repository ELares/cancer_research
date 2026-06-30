{pkgs ? import <nixpkgs> {}}:
let
  python = pkgs.python312;
  tex = pkgs.texlive.combine {
    inherit (pkgs.texlive)
      scheme-medium
      latexmk
      collection-fontsrecommended
      collection-latexextra;
  };
in
pkgs.mkShell {
  packages = with pkgs; [
    python
    uv
    cargo
    rustc
    rustfmt
    clippy
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

    echo "Entered cancer_research nix-shell fallback"
    echo "Prefer \`nix develop\` for the exact Rust toolchain pin from simulations/rust-toolchain.toml."
    echo "Run: uv sync --frozen"
  '';
}
