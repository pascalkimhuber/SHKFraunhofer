(TeX-add-style-hook
 "lennard_jones"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "scrartcl"
    "scrartcl10"
    "mathtools"
    "amsmath"
    "amssymb")
   (TeX-add-symbols
    "boldp"
    "boldq")
   (LaTeX-add-labels
    "sec:hess-lenn-jones"
    "eq:1"
    "eq:2"
    "sec:part-deriv-r"
    "sec:one-dimens-deriv"
    "sec:lennard-jones-forces")))

