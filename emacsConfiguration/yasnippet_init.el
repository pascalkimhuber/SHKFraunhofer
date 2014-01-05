;;; Load and enable yasnippet mode. 
(require 'yasnippet)
(yas-global-mode t) 

;;; Add snippet directories to snippet tables. 
(yas/load-directory "~/.emacs.d/elpa/yasnippet-20131021.928/snippets")


;;; Load snippet tables. 
(yas-reload-all)

;;; Enable drowpdown lists. 
(require 'dropdown-list)

;;; Set prompting method. 
(setq yas/prompt-functions '(yas/dropdown-prompt yas/ido-prompt))
