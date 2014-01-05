;;; Init file for spell checking configuration (flyspell and flymake). 

;; Bind 'C-c f' to flyspell-mode. 
(global-set-key (kbd "C-c f") 'flyspell-mode)

;; Change dictionary automaticaly. This uses the package auto-dictionary. 
(require 'auto-dictionary)
(add-hook 'flyspell-mode-hook (lambda () (auto-dictionary-mode 1)))

;; TODO Add configuration for flymake-mode. 
