;;=== Package manager =========================================================

;;; Load package sources.

;;; ELPA is the Emacs Lisp Package Archive.
;;; MELPA is the Milkypostman's Emacs Lisp Package Archive.
(require 'package)
(dolist (source '(("marmalade" . "http://marmalade-repo.org/packages/")
		  ("elpa" . "http://tromey.com/elpa/")
		  ("melpa" . "http://melpa.milkbox.net/packages/")
		  ))
  (add-to-list 'package-archives source t))
(package-initialize)
