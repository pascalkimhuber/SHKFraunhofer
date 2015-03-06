(TeX-add-style-hook "makefile"
 (lambda ()
    (LaTeX-add-labels
     "sec:writing-makefiles"
     "sec:targ-rules-depend"
     "sec:bemerkung"
     "sec:beispiel"
     "sec:das-defaulttarget"
     "sec:beispiel-1"
     "sec:pattern-regeln"
     "sec:beispiel-2"
     "sec:variablen-makefiles"
     "sec:beispiel-3"
     "sec:kommentare-makefiles"
     "sec:phony-targets"
     "sec:beispiel-4"
     "sec:pattern-substitution"
     "sec:beispiel-5"
     "sec:abhang-als-targ"
     "sec:beispiel-6"
     "sec:rekursives-make")
    (TeX-run-style-hooks
     "color"
     "listings"
     "inputenc"
     "utf8"
     "latex2e"
     "scrartcl10"
     "scrartcl"
     "")))

