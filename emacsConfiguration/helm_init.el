;;; TODO Learn how to use helm. Add usefull configuration. 

;;; Enable helm-mode and load required helm-libraries. 
(require 'helm-config)
(require 'helm-buffers)
(require 'helm-eval)
(require 'helm-imenu)
(require 'helm-net)
(require 'helm-files)

;;; Enable recursive-minibuffers.
(setq enable-recursive-minibuffers t)

;;; Redefine helm-mini. 
; Enable follow-mode for multi-occur, buffer-list and browse-code.
(eval-after-load "helm-regexp"
  '(helm-attrset 'follow 1 helm-source-buffers-list))

; Redefine helm-mini and bind it to C-c h (for helm). 
(defun helm-mini ()
  "Preconfigured `helm' lightweight version \(buffer -> recentf\)."
  (interactive)
  (let ((buffers (list (current-buffer))))
    (helm-other-buffer '(helm-c-source-buffers-list
			 helm-c-source-google-suggest
			 helm-c-source-calculation-result
			 helm-c-source-find-files
			 helm-c-source-recentf     
			 helm-c-source-info-pages
			 helm-c-source-buffer-not-found
			 helm-c-source-imenu)
		       "*helm mini*")))
(global-set-key (kbd "C-c h") 'helm-mini)
