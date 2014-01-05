;;; Enable smart-mode-line. 
(require 'smart-mode-line)
(if after-init-time (sml/setup)
  (add-hook 'after-init-hook 'sml/setup))

;;; Show column number in mode-line. 
(setq column-number-mode t)

;;; Customize colors.
(setq
 sml/active-background-color "gray60"
 sml/active-foreground-color "gray40")

;;; Customize faces.
(set-face-attribute 'sml/global nil
 		    :foreground "gray40")
(set-face-attribute 'sml/col-number nil
		    :inherit 'sml/global 
		    :foreground "sea green")
(set-face-attribute 'sml/filename nil
		    :inherit 'sml/global 
		    :foreground "gray0" 
		    :weight 'bold)
(set-face-attribute 'sml/line-number nil
		    :inherit 'sml/global 
		    :foreground "blue2" 
		    :weight 'bold)
(set-face-attribute 'sml/modes nil
		    :inherit 'sml/global 
		    :foreground "dark red")
(set-face-attribute 'sml/position-percentage nil
		    :inherit 'sml/filename 
		    :foreground "DarkGoldenrod4" 
		    :weight 'normal)
