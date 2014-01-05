;;; Enable and customize auto-complete. 
; The manual can be found at http://cx4a.org/software/auto-complete/manual.html#auto-complete_command.

;;; Add auto-complete directory to load-path.
(add-to-list 'load-path "~/.emacs.d/elpa/auto-complete")

;;; Basic configuration.
(require 'auto-complete-config)
(add-to-list 'ac-dictionary-directories "~/.emacs.d/elpa/auto-complete/dict")
(ac-config-default)

;;; Bind auto-complete command to "M-TAB".
(define-key ac-mode-map (kbd "M-TAB") 'auto-complete)

;;; Use ac flyspell workaround. 
(ac-flyspell-workaround)

;;; Select auto-completion candidates with 'C-n' and 'C-p'. 
(setq ac-use-menu-map t)
(define-key ac-menu-map "\C-n" 'ac-next)
(define-key ac-menu-map "\C-p" 'ac-previous)

;;; Set sources in cc-mode (using ac-clang-complete). (See manual for details.)
(require 'auto-complete-clang-async)

(defun ac-cc-mode-setup ()
  (setq ac-clang-complete-executable "~/.emacs.d/clang-complete")
  (setq ac-sources '(ac-source-clang-async))
  (ac-clang-launch-completion-process)
)

(defun my-ac-config ()
  (add-hook 'c-mode-common-hook 'ac-cc-mode-setup)
  (add-hook 'auto-complete-mode-hook 'ac-common-setup)
  (global-auto-complete-mode t))

(my-ac-config)

(add-to-list 'load-path "~/.emacs.d/auto-complete-clang/")
(require 'auto-complete-clang)

