J/A+A/704/A344        eRASS1 counterparts and classification    (Salvato+, 2025)
================================================================================
Counterpart identification and classification for eRASS1 and characterisation of
the AGN content.
    Salvato M., Wolf J., Dwelly T., Starck H., Buchner J., Shirley R.,
    Merloni A., Georgakakis A., Balzer F., Brusa M., Rau A., Freund S., Lang D.,
    Liu T., Lamer G., Schwope A., Roster W., Waddell S., Scialpi M., Igo Z.,
    Kluge M., Mannucci F., Tiwari S., Homan D., Krumpe M., Zenteno A.,
    Hernandez-Lang D., Comparat J., Fabricius M., Snigula J., Schlegel D.,
    Weaver B.A., Zhou R., Dey A., Valdes F., Myers A., Juneau S., Winkler H.,
    Marquez I., di Mille F., Ciroi S., Schramm M., Buckley D.A.H., Brink J.,
    Gromadzki M., Robrade J., Nandra K.
    <Astron. Astrophys. 704, A344 (2025)>
    =2025A&A...704A.344S        (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: QSOs ; Galaxies ; X-ray sources ; Optical ; Infrared
Keywords: methods: data analysis - methods: statistical - catalogs - surveys -
          galaxies: active - X-rays: general

Abstract:
    Accurately accounting for the AGN phase in galaxy evolution requires a
    large, clean AGN sample. This is now possible with
    SRG/eROSITA, which completed its first all-sky X-ray survey (eRASS1)
    on June 12, 2020. The public Data Release 1 (DR1, Jan 31, 2024)
    includes 930,203 sources from the Western Galactic Hemisphere.

    The data enable the selection of a large AGN sample and the discovery
    of rare sources. However, scientific return depends on accurate
    characterisation of the X-ray emitters, requiring high-quality
    multiwavelength data. This paper presents the identification and
    classification of optical and infrared counterparts to eRASS1 sources.

    Counterparts to eRASS1 X-ray point sources were identified using Gaia
    DR3, CatWISE2020, and Legacy Survey DR10 (LS10) with the Bayesian NWAY
    algorithm and trained priors. Sources were classified as Galactic or
    extragalactic via a Machine Learning model combining optical/IR and
    X-ray properties, trained on a reference sample. For extragalactic
    LS10 sources, photometric redshifts were computed using Circlez.

    Within the LS10 footprint, all 656,614 eROSITA/DR1 sources have at
    least one possible optical counterpart; ~570,000 are extragalactic and
    likely AGN. Half are new detections compared to AllWISE, Gaia, and
    Quaia AGN catalogues. Gaia and CatWISE2020 counterparts are less
    reliable, due to the survey's shallowness and the limited amount of
    features available to assess the probability of being an X-ray
    emitter. In the Galactic Plane, where the overdensity of stellar
    sources also increases the chance of associations, using conservative
    reliability cuts, we identify approximately 18000 Gaia and 55000
    CatWISE2020 extragalactic sources.

    We release three high-quality counterpart catalogues - plus the
    training and validation sets - as a benchmark for the field. These
    datasets have many applications, but in particular empower researchers
    to build AGN samples tailored for completeness and purity,
    accelerating the hunt for the Universe's most energetic engines.

Description:
    In this paper, we describe the construction of the three catalogues
    listing the most reliable counterparts to the point sources
    detected by eROSITA (PSF Half Energy Width of ~30" and corresponding
    1{sigma} median positional uncertainty of ~4.5" and a uniform flux
    limit at 50% completeness of F0.5-2keV>5x10^-14^erg/s/cm^2^) in the
    first all-sky pass over the Western hemisphere (eROSITA_DE), using
    LS10, Gaia DR3, and CW2020, respectively.

    All catalogues constructed in this paper are released to the
    public (see Section 8) via Vizier and the eROSITA web pages at
    https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/
    (files 4-12).

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
dr1ls10.dat     2625   966194   Identification and characterization of the
                                 counterparts to the X-ray point sources in the
                                 eROSITA/DR1 main catalogue using the
                                 Legacy Survey DR10
dr1ge3.dat      1671   959261   Identification and characterization of the
                                 counterparts to the X-ray point sources in the
                                 eROSITA/DR1 main catalogue using the
                                 Gaia DR3 survey
dr1cw2020.dat   1456  1001527   Identification and characterization of the
                                 counterparts to the X-ray point sources in the
                                 eROSITA/DR1 main catalogue using the
                                 CatWISE2020 survey
csc2ls10.dat     220   640567   Legacy Survey DR10 counterparts to Chandra CSC2
                                 and field population used for training
4xmmls10.dat     210  4140639   Legacy Survey DR10 counterparts to 4XMM and
                                 field population used for training
csc2ge3.dat      212   143315   Gaia EDR3 counterparts to Chandra CSC2 and
                                 field population used for training
4xmmge3.dat      202   440544   Gaia EDR3 counterparts to 4XMM and field
                                 population used for training
csc2wise.dat     211   201404   CatWISE2020 counterparts to Chandra CSC2 and
                                 field population used for training
4xmmcw20.dat     201  1328859   CatWISE2020 counterparts to 4XMM and field
                                 population used for training
--------------------------------------------------------------------------------

See also:
 J/A+A/682/A34 : SRG/eROSITA all-sky survey catalogs (eRASS1) (Merloni+, 2024)

 IX/69  : XMM-Newton Serendipitous Source Catalogue 4XMM-DR13 (Webb+, 2023)
 IX/70  : Chandra Source Catalog Release 2 (CSC 2.1) (Evans+, 2024)
 I/355  : Gaia DR3 Part 1. Main source (Gaia Collaboration, 2022)
 II/365 : The CatWISE2020 catalog (updated version 28-Jan-2021) (Marocco+, 2021)

 https://www.legacysurvey.org/dr10 : DESI LS10 Home Page

Byte-by-byte Description of file: dr1ls10.dat
--------------------------------------------------------------------------------
   Bytes   Format Units        Label        Explanations
--------------------------------------------------------------------------------
    1-  23  A23   ---          IAUName      eROSITA String containing the
                                             official IAU name of the source
                                             (IAUNAME)
   25-  56  A32   ---          DetUId       eROSITA String unique detection
                                             (DETUID)
   58-  63  I6    ---          SkyTile      eROSITA Sky tile ID (SKYTILE)
   65-  68  I4    ---          IdSrc        eROSITA Source ID in each sky tile
                                             (ID_SRC) (1)
   70-  81  I12   ---          UId          eROSITA Integer unique detection ID
                                             (UID) (2)
   83-  95  I13   ---          UIdHard      eROSITA Hard catalog UID of the
                                             source with a strong association,
                                             or-UID if the association is weak
                                             (UID_Hard)
   97- 100  I4    ---          IdCluster    eROSITA Group ID of simultaneously
                                             fitted sources (ID_CLUSTER)
  102- 116 F15.11 deg          RAdeg        eROSITA Right ascension (ICRS),
                                             corrected (RA)
  118- 132 F15.11 deg          DEdeg        eROSITA Declination (ICRS),
                                             corrected (DEC_1)
  134- 148 F15.11 deg          RARdeg       eROSITA Right ascension (ICRS),
                                             uncorrected (RA_RAW)
  150- 164 F15.11 deg          DERdeg       eROSITA Declination (ICRS),
                                             uncorrected (DEC_RAW)
  166- 175  F10.6 arcsec     e_RAdeg        ?=- eROSITA 1-sigma lower error on
                                             RA (RA_LOWERR)
  177- 186  F10.6 arcsec     E_RAdeg        ?=- eROSITA 1-sigma upper error on
                                             RA (RA_UPERR)
  188- 197  F10.6 arcsec     e_DEdeg        ?=- eROSITA 1-sigma lower error on
                                             DEC (DEC_LOWERR)
  199- 208  F10.6 arcsec     E_DEdeg        ?=- eROSITA 1-sigma upper error on
                                             DEC (DEC_UPERR)
  210- 219  F10.6 arcsec       RADEErr      eROSITA Combined positional error,
                                             raw output from PSF fitting
                                             (RADEC_ERR)
  221- 230  F10.6 ---          posErr       eROSITA 1-sigma positional
                                            uncertainty (POS_ERR)
  232- 246 F15.11 deg          GLON         Galactic longitude (LII)
  248- 262 F15.11 deg          GLAT         Galactic latitude (BII)
  264- 278 F15.11 deg          ELON         Ecliptic longitude (ELON)
  280- 294 F15.11 deg          ELAT         Ecliptic latitude (ELAT)
  296- 302  F7.1  d            MJD          Modified Julian Date of the
                                             observation of the source nearest
                                             to the optical axis of the
                                             telescope (MJD)
  304- 310  F7.1  d          b_MJD          Minimum Modified Julian Date of
                                             observations used to construct the
                                             model of this object (MJD_MIN)
  312- 318  F7.1  d          B_MJD          Maximum Modified Julian Date of
                                             observations used to construct the
                                             model of this object (MJD_MAX)
  320- 327  F8.6  arcsec       Ext          [0] Source extent (EXT)
  329- 336  F8.6  arcsec     s_Ext          [0]?=- Extent error (EXT_ERR)
  338- 345  F8.6  arcsec     e_Ext          [0]?=- 1-sigma lower error on EXT
                                             (EXT_LOWERR)
  347- 354  F8.6  arcsec     E_Ext          [0]?=- 1-sigma upper error on EXT
                                             (EXT_UPERR)
  356- 363  F8.6  ---          ExtLike      [0] Extent likelihood (EXT_LIKE)
  365- 376  E12.6 ---          DetLike0     Detection likelihood measured by
                                             forced PSF-fitting at 0.2-2.3keV
                                             (DET_LIKE_0_1)
  378- 389  E12.6 ct           MLcts1       Source net counts measured from
                                             count rate at 0.2-2.3keV (ML_CTS_1)
  391- 402  E12.6 ct         s_MLcts1       ?=- 1 sigma counts error at
                                             0.2-2.3keV (ML_CTS_ERR_1)
  404- 415  E12.6 ct         e_MLcts1       ?=- 1 sigma lower error of counts at
                                             0.2-2.3keV (ML_CTS_LOWERR_1)
  417- 428  E12.6 ct         E_MLcts1       ?=- 1 sigma upper error of counts at
                                             0.2-2.3keV (ML_CTS_UPERR_1)
  430- 439  F10.6 ct/s         MLRate1      ?=- Source count rate measured by
                                             forced PSF-fitting; 0.2-2.3keV
                                             (ML_RATE_1)
  441- 449  F9.6  ct/s       s_MLRate1      ?=- 1 sigma count rate erro