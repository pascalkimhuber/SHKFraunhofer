;;; Init file for multiple-cursor configuration. A full list of keys can be found at https://github.com/magnars/multiple-cursors.el. 

;; Bind C-c c to 'mc/edit-lines'. 
(global-set-key (kbd "C-c c") 'mc/edit-lines)

;; Bind C-> to 'mc/mark-next-like-this'. 
(global-set-key (kbd "C->") 'mc/mark-next-like-this)

;; Bind C-< to 'mc/mark-previous-like-this'. 
(global-set-key (kbd "C-<") 'mc/mark-previous-like-this)

;; Vind C-. to 'mc/mark-all-like-this'. 
(global-set-key (kbd "C-.") 'mc/mark-all-like-this)
