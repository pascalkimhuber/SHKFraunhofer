(TeX-add-style-hook "Tremolo_Manual"
 (lambda ()
    (LaTeX-add-index-entries
     "TODO: #2"
     "TODO: #1"
     "FIXME: #2"
     "FIXME: #1")
    (LaTeX-add-bibliographies
     "./bib/Tremolo_Manual.bib")
    (TeX-add-symbols
     '("fixme" ["argument"] 1)
     '("todo" ["argument"] 1)
     '("draft" 1)
     '("fehler" 1)
     '("notiz" 1)
     "angstrom"
     "astronomicalunit"
     "siday"
     "siyear"
     "laplace")
    (TeX-run-style-hooks
     "makeidx"
     "draftwatermark"
     "ifthen"
     "pstricks"
     "booktabs"
     "siunitx"
     "listings"
     "color"
     "framed"
     "hyperref"
     "amssymb"
     "amsmath"
     "babel"
     "english"
     "graphicx"
     "pdftex"
     "fontenc"
     "T1"
     "inputenc"
     "utf8x"
     "latex2e"
     "bk11"
     "book"
     "a4paper"
     "11pt"
     "./Content/Overview"
     "./Content/ConfigurationSetup"
     "./Content/FirstSteps"
     "./Content/RunningTremolo"
     "./Content/DomainAndBoundaries"
     "./Content/EnsemblesAndThermostats"
     "./Content/Optimization"
     "./Content/Longrange"
     "./Content/Parallelization"
     "./Content/MeasurementAndOutput"
     "./Content/Potentials"
     "./Content/External"
     "./Content/InputFiles"
     "./Content/Tutorial")))

