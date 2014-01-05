;;; Enable and customize org-mode. 
; org-mode supports:
; - sparse trees
; - plain lists
; - footnotes
; - tables
; - hyperlinks (external and internal)
; - TODO items (with process logging, priorities, subtasks, checkboxes)
; - tags
; - properties 
; - timestamps (deadlines etc.)
; - capture (for capturing new ideas), archiving
; - agenda views (sort of overview)
; - export of org-files
; - latex support
; - publishing
; - source code support
; - mobileOrg

;;; Enable org-mode for .org, .org_archive and .txt files by default.
(add-to-list 'auto-mode-alist '("\\.\\(org\\|org_archive\\)$" . org-mode))

;;; Set global keys for the most important org commands. 
(global-set-key "\C-cl" 'org-store-link)
(global-set-key "\C-cc" 'org-capture)
(global-set-key "\C-ca" 'org-agenda)
(global-set-key "\C-cb" 'org-iswitchb)

;;; TODO Customize sequential states of todos (see 5.2 in org-guide). 

;;; TODO Customize capture templates (see 9.1.3 in org-guide).

