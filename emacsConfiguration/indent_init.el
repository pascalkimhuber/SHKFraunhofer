;;; Customize indentation. 

;;; Enable tabs for indentation (this is the default in emacs). 
(setq-default indent-tabs-mode t)

;;; Set tab-width to 4 and equal to indentation offset in CC Mode and Perl. 
(setq-default tab-width 4)

;;; Enable smart-tabs-mode package in some major modes. 
(smart-tabs-insinuate 'c 'c++ 'javascript 'java 'python 'ruby)

;;; Enable k&r style for c-like modes. (See http://www.emacswiki.org/emacs/IndentingC for details). 
(setq c-default-style "k&r"
      c-basic-offset 4)
