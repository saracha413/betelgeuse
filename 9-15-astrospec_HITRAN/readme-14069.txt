Basic
-----

The fields output for this profile are listed in order below with format strings, units and description the following information: :

nu
--
Data type: float
Units: cm-1
Description: Transition wavenumber

sw
--
Data type: float
Units: cm-1/(molec.cm-2)
Description: Line intensity, multiplied by isotopologue abundance, at T = 296 K

molec_id
--------
Data type: int
Units: [dimensionless]
Description: The HITRAN integer ID for this molecule in all its isotopologue forms

local_iso_id
------------
Data type: int
Units: [dimensionless]
Description: Integer ID of a particular Isotopologue, unique only to a given molecule, in order or abundance (1 = most abundant)

nu-err
------
Data type: int
Units: [dimensionless]
Description: HITRAN uncertainty code for nu

sw-err
------
Data type: int
Units: [dimensionless]
Description: HITRAN uncertainty code for sw

nu-ref
------
Data type: int
Units: [dimensionless]
Description: Source (reference) ID for nu

sw-ref
------
Data type: int
Units: [dimensionless]
Description: Source (reference) ID for sw
