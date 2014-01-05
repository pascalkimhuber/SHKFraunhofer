;;; TODO: Install email package for reading mails in emacs.

;;; Sending Mail with Emacs [Gmail] using SMTP.
(setq 
 user-full-name "Pascal Huber"
 user-mail-address "pascalkimhuber@gmail.com")

(setq 
 smtpmail-stream-type 'ssl
 send-mail-function 'smtpmail-send-it
 smtpmail-smtp-server "smtp.gmail.com"
 smtpmail-smtp-service 465)

(require 'smtpmail)



