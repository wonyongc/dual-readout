{
  description = "dual-readout";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/22.11";
    utils.url = "github:numtide/flake-utils";
    hepnix.url = "github:wonyongc/hepnix";
    hepnix.inputs.nixpkgs.follows = "nixpkgs";
  };

  outputs = { self, nixpkgs, utils, hepnix, ... }@inputs:

    inputs.utils.lib.eachSystem [ "x86_64-linux" "aarch64-darwin" ]

    (system: let
      pkgname = "dual-readout";

#      hepcore = hepnix.outputs.hepcore;

      pkgs = import nixpkgs {
        inherit system;
        overlays = [];
        config.allowUnfree = true;
      };

      in {

        packages.default = pkgs.stdenv.mkDerivation rec {
          pname = pkgname;
          version = "flake-test";

          src = self;

          buildInputs = inputs.hepnix.hepcore.core nixpkgs ++ [
            hepnix.heppkgs.podio
            hepnix.heppkgs.edm4hep
            hepnix.heppkgs.VecCore
            hepnix.heppkgs.VecGeom
            hepnix.heppkgs.vdt
            hepnix.heppkgs.SIO
            hepnix.heppkgs.LCIO
            hepnix.heppkgs.geant4
            hepnix.heppkgs.dd4hep
            hepnix.heppkgs.evtgen
            hepnix.heppkgs.gaudi
            hepnix.heppkgs.K4FWCore
            hepnix.heppkgs.k4SimGeant4
            hepnix.heppkgs.k4Gen ];

          nativeBuildInputs = hepnix.hepcore.wrappers nixpkgs;
        };
      }
    );
}