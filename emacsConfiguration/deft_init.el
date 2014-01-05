;;; Enable the deft file browser (see http://jblevins.org/projects/deft/). 
; Notes are saved in ~/.deft by default (this can be changed by: 
;(setq deft-directory "/path/to/directory")). 
(require 'deft)
(setq deft-directory "~/.org")

;;; Use org-files for deft-notes. 
(setq deft-extension "org")
(setq deft-text-mode 'org-mode)

;;; Bind C-c d (for deft) to deft-mode. 
(global-set-key (kbd "C-c d") 'deft)




