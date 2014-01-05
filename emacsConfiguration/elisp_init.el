;;; Configuration for Emacs-Lisp mode. 

;;; Set auto-complete configuration (add some sources). 
(defun ac-elisp-mode-setup ()       
  (setq ac-sources
	(append '(ac-source-filename
		  ac-source-files-in-current-dir)
		ac-sources)))
(add-hook 'emacs-lisp-mode-hook 'ac-elisp-mode-setup)




