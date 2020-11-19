# isotopic-enrichment-calculator

Isotopic composition of the constituent elements dictates the shape of isotope patterns observed in mass spectrometry. As an example, a molecule with elemental formula C_{c}H_{h}N_{n}O_{o} will contain (1 + c)(1 + h)(1 + n)(1 + o)(2 + o)/2 distinct combinations of ions (called isotopologues). 

The abundance (proportion) of all of these isotopologues is governed in accordance with the basic rules of probability. As an example, the lightest of all the ions will have an approximate mass of c*m(12C) + h*m(1H) + n*m(14N) + o*m(16O) and will have an abundance given by the followong expression: x(12C)^c * x(1H)^h * x(14N)^n * x(16O)^o, where x(12C) is the abundance (proportion) of carbon-12 among all carbon atoms in the molecule. Thus, the knowledge of the molecular formula and the isotopic composition of all makeup elements enables us to establish the expected (‘theoretical’) isotope patterns of molecules.

This calculator performs the reverse calculation: knowing the identity of a compound (i.e. knowing its molecular formula) and by having an observed isotopic pattern of it, one can calculate the isotopic composition of one of its constituent atoms that best matches the observed isotope pattern.

This calculator was developed to study the incorporation of nitrogen-15 in marine biotoxins and was described in a peer-reviewed publication available at https://doi.org/10.3390/md17110643.
