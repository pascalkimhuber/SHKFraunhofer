;;; TODO Add all latex functionalities to emacs. This includes:
;;; AUCTeX (advanced LaTeX mode)
;;; RefTeX (make citations, refereces, table of contents, labels, and more)

;;; Add latex executables to exec path. 
(setenv "PATH" (concat (getenv "PATH") ":/usr/texbin"))
    (setq exec-path (append exec-path '("/usr/texbin")))

;;; Add pdf2dsc (ghostscript) executables to exec path, used by preview-latex.
(setenv "PATH" (concat (getenv "PATH") ":/opt/local/bin"))
    (setq exec-path (append exec-path '("/opt/local/bin")))

;;; Enable document parsing.
(setq TeX-auto-save t)
(setq TeX-parse-self t)

;;; Enable multi-file document structure by default.
(setq-default TeX-master nil)

;;; Enable some minor modes. 
(add-hook 'LaTeX-mode-hook 'visual-line-mode)    ;; enable line wrapping
(add-hook 'LaTeX-mode-hook 'flyspell-mode)       ;; enable spell checking
(add-hook 'LaTeX-mode-hook 'LaTeX-math-mode)     ;; enable math mode

;;; Compile documents to PDF by default. 
(setq TeX-PDF-mode t)

;;; Enable forward and inverse search with Viewer.
;;; In order to use inverse search press Cmd+Shift and click on line in Skim. 
(setq TeX-source-correlate-mode t)
(setq TeX-source-correlate-start-server (quote ask))
(setq TeX-source-correlate-start-server t)

;;; Use Skim as PDF Viewer.
(setq TeX-view-program-selection
      '((output-dvi "DVI Viewer")
        (output-pdf "PDF Viewer")
        (output-html "HTML Viewer")))
(setq TeX-view-program-list
      '(("DVI Viewer" "open %o")
        ("PDF Viewer" "/Applications/Skim.app/Contents/SharedSupport/displayline -g %n %o %b")
        ("HTML Viewer" "open %o")))

;;; Enable RefTeX mode for references, citations etc. 
(add-hook 'LaTeX-mode-hook 'turn-on-reftex)
(setq reftex-plug-into-AUCTeX t)

;;; Enable preview-latex by setting Gs options (important: -dNOSAFER).
(setq preview-gs-options '("-q" "-dNOSAFER" "-dNOPAUSE" "-DNOPLATFONTS" "-dPrinted" "-dTextAlphaBits=4" "-dGraphicsAlphaBits=4")) 

;;; Enable yasnippets. 
(add-hook 'LaTeX-mode-hook 
	  '(lambda ()
	     (yas-minor-mode)))

;;; Enable auto-completion in LaTeX-mode. 
(require 'ac-math)
; Make auto-complete aware of 'LaTeX-mode'.
(add-to-list 'ac-modes 'LaTeX-mode) 
; Add ac-sources to default ac-sources.
(defun ac-latex-mode-setup ()       
  (setq ac-sources
	(append '(ac-source-math-unicode 
		  ac-source-math-latex 
		  ac-source-latex-commands
		  ac-source-yasnippet)
		ac-sources)))
(add-hook 'LaTeX-mode-hook 'ac-latex-mode-setup)
; Activate auto-completion in LaTeX-mode. 
(add-hook 'LaTeX-mode-hook 'auto-complete-mode)

