# CANDY (Chemistry ANd DYnamics)

CANDY contains two submodules, a modified astrochem (Maret & Bergin 2015) directory and the chemdiff directory, where the python wrapper is contained.

## Astrochem

The astrochem folder is copied from the (astrochem)[https://github.com/smaret/astrochem] github page, complete with documentation and installation instructions. Changes have been made within the `src` folder, and Astrochem can be installed as usual. If you already have astrochem installed on your machine, you can simply copy the `src` folder from here into your astrochem directory (replacing the default astrochem/src folder), then reinstall as usual.

Chemical network files (`.chm`) should follow the same format as astrochem but with added chemical reactions for:
	|description | reacno|
	|------------:|:-------|
	| HD formation on grains | 100 |
	| Shielded dissociation of H2 | 14 |
	| Shielded dissociation of HD | 114 |
	| Shielded photo-ionization of C | 15 |
	| Shielded dissociation of C16O | 16 |
	| Shielded dissociation of C17O | 17 |
	| Shielded dissociation of C18O | 18 |
	| Hydrogenation | 25 |
	| Cosmic-ray desorption of CO | 42 |
	| Cosmic-ray reactions involving He | 43 |
	| secondary xray ionization of ... H | 60 |
	| ... H2 | 61 |
	| ... other molecules | 62 |
	| creation of excited H2 | 90 |
	| de-excitation of H2 | 91 | 
	| reactions with excited H2 | 92 |