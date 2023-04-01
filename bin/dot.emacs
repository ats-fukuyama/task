;
(autoload 'f90-mode "f90"
 "Major mode for editing Fortran 90 code in free format." t)
 (setq auto-mode-alist (append auto-mode-alist 
                       (list '("\\.f90$" . f90-mode))))
;
(set-default-coding-systems 'utf-8)
(set-buffer-file-coding-system 'utf-8)
;(set-terminal-coding-system 'euc-jp)
(set-terminal-coding-system 'utf-8)
;(set-keyboard-coding-system 'euc-jp)
(set-keyboard-coding-system 'utf-8)
;; inhibit startup message
(setq inhibit-startup-message t)
;; no toolbar
;; (tool-bar-mode -1)
(menu-bar-mode -1)
;; no scratch buffer message
(setq initial-scratch-message "")
;;
;; 
(setq fortran-comment-indent-style nil)
(setq blink-matching-delay 0.1)
(put 'eval-expression 'disabled nil)
(put 'upcase-region 'disabled nil)

;(load "~/lib/emacs/ws8.el")
(custom-set-variables
 ;; custom-set-variables was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(column-number-mode t)
 '(tool-bar-mode nil nil (tool-bar))
 '(warning-suppress-types '((comp))))
(custom-set-faces
 ;; custom-set-faces was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 )

(put 'downcase-region 'disabled nil)

(add-to-list 'auto-mode-alist '("\\.F90\\'" . f90-mode))

(setq-default ispell-program-name "aspell")
(with-eval-after-load "ispell"
  (setq ispell-local-dictionary "en_US")
  (add-to-list 'ispell-skip-region-alist '("[^\000-\377]+")))
(add-hook 'latex-mode-hook 'flyspell-mode)

(setq-default truncate-lines nil)
(setq-default trancate-partial-width-windows nil)
(setq truncate-lines nil)
(setq trancate-partial-width-windows nil)
;;(global-set-key (kbd "C-x t") 'toggle-truncate-lines)
