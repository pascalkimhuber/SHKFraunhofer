;;; Enable ido-mode for fast buffer and file management.
(setq ido-everywhere t)
(ido-mode 1)

;;; Enable flex ("fuzzy") matching. 
(setq ido-enable-flex-matching t)

;;; Disable prompts for non-existent buffers.
(setq ido-create-new-buffer 'always)
(setq confirm-nonexistent-file-or-buffer nil)
