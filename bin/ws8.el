;;
;; Wordstar emulation setting  V6.00  2004/04/23
;;
;; Author: A. Fukuayama <fukuyama@nucleng.kyoto-u.ac.jp>
;;
;; Part of the code is borrowed from 
;;    ws-mode.el --- WordStar emulation mode for GNU Emacs
;;    Copyright (C) 1991 Free Software Foundation, Inc.
;;    Author: Juergen Nickelsen <nickel@cs.tu-berlin.de>
;;    Version: 0.7
;;

;; interrupt character has been changed to ^L

(set-input-mode nil nil 0 ?\^l)
;; (electric-indent-mode 0)

;; Key-in character table

;(let ((the-table (make-string 128 0)))
;  (let ((i 0))
;    (while (< i 128)
;      (aset the-table i i)
;      (setq i (1+ i))))
;; Swap ^L and ^g
;  (aset the-table ?\^l ?\^g)
;  (aset the-table ?\^g ?\^l)
;; Swap ^C and ^V
;  (aset the-table ?\^c ?\^v)
;  (aset the-table ?\^v ?\^c)
;; Swap ^X and ^K
;  (aset the-table ?\^x ?\^k)
;  (aset the-table ?\^k ?\^x)
;; transrate ^^ as ^s and ^_ as ^q for XON/XOFF
;  (aset the-table ?\^^ ?\^s)
;  (aset the-table ?\^_ ?\^q)
;  (setq keyboard-translate-table the-table))

(keyboard-translate ?\^l ?\^g)
(keyboard-translate ?\^g ?\^l)
(keyboard-translate ?\^c ?\^v)
(keyboard-translate ?\^v ?\^c)
(keyboard-translate ?\^k ?\^x)
(keyboard-translate ?\^x ?\^k)
(keyboard-translate ?\^^ ?\^s)
(keyboard-translate ?\^_ ?\^q)

(setq scroll-step 1)
;;(setq fortran-line-number-indent 5)
;;(setq fortran-continuation-string "&")
;;(setq fortran-check-all-num-for-matching-do t)

(define-key isearch-mode-map "\C-f" 'isearch-repeat-forward)
(define-key isearch-mode-map "\C-a" 'isearch-repeat-backward)
(define-key isearch-mode-map "\C-s" 'isearch-other-control-char)
(define-key isearch-mode-map "\C-r" 'isearch-other-control-char)
(define-key isearch-mode-map "\C-h" 'isearch-delete-char)
(define-key isearch-mode-map "\C-p" 'isearch-quote-char)
(define-key isearch-mode-map "\C-q" 'isearch-other-control-char)

(define-key minibuffer-local-isearch-map "\C-f" 'isearch-forward-exit-minibuffer)
(define-key minibuffer-local-isearch-map "\C-a" 'isearch-reverse-exit-minibuffer)
(define-key minibuffer-local-isearch-map "\C-s" nil)
(define-key minibuffer-local-isearch-map "\C-r" nil)

(define-key edit-abbrevs-map "\C-k" nil)
(define-key edit-abbrevs-map "\C-x\C-r" 'edit-abbrevs-redefine)

;;; global-map

(define-key global-map "\C-a" 'backward-word)
(define-key global-map "\C-b" 'fill-paragraph)
(define-key global-map "\C-v" 'scroll-up)                        ;\C-c
(define-key global-map "\C-d" 'forward-char)
(define-key global-map "\C-e" 'previous-line)
(define-key global-map "\C-f" 'forward-word)
(define-key global-map "\C-l" 'delete-char)                      ;\C-g
(define-key global-map "\C-h" 'backward-delete-char-untabify)
;;(define-key global-map "\C-i" 'indent-for-tab-command)
;;(define-key global-map "\C-j" 'help-for-help)
(define-key global-map "\C-x" 'Control-X-prefix)
(define-key global-map "\C-g" 'keyboard-quit) ; \C-l
;;                      \C-m   newline
(define-key global-map "\C-n" 'universal-argument)
(define-key global-map "\C-o" 'ctrl-O-prefix)
(define-key global-map "\C-p" 'quoted-insert)
(define-key global-map "\C-q" 'ctrl-Q-prefix)
(define-key global-map "\C-r" 'scroll-down)
(define-key global-map "\C-s" 'backward-char)
(define-key global-map "\C-t" 'kill-word)
(define-key global-map "\C-u" 'undo)
(define-key global-map "\C-c" 'mode-specific-command-prefix)     ;\C-v
(define-key global-map "\C-w" 'scroll-down-line)
(define-key global-map "\C-k" 'next-line)
(define-key global-map "\C-y" 'kill-complete-line)
(define-key global-map "\C-z" 'scroll-up-line)
;;                      \C-[   Prefix Command
;; (define-key global-map "\C-\\" nil) ;; reserved for egg
;; (define-key global-map "\C-]" nil) ;; reserved telnet escape
;; (define-key global-map "\C-^"  nil) ;; \C-s
;; (define-key global-map "\C-_"  nil) ;; \C-q
;; (define-key global-map "\C-?" 'backward-delete-char-untabify)

;; C-k-map
  
(define-key ctl-x-map " " ())
(define-key ctl-x-map "0" 'ws-set-marker-0)
(define-key ctl-x-map "1" 'ws-set-marker-1)
(define-key ctl-x-map "2" 'ws-set-marker-2)
(define-key ctl-x-map "3" 'ws-set-marker-3)
(define-key ctl-x-map "4" 'ws-set-marker-4)
(define-key ctl-x-map "5" 'ws-set-marker-5)
(define-key ctl-x-map "6" 'ws-set-marker-6)
(define-key ctl-x-map "7" 'ws-set-marker-7)
(define-key ctl-x-map "8" 'ws-set-marker-8)
(define-key ctl-x-map "9" 'ws-set-marker-9)
(define-key ctl-x-map "a" 'ws-append-to-buffer)
(define-key ctl-x-map "\C-a" 'ws-append-to-file)
(define-key ctl-x-map "\C-b" 'ws-set-begining-of-block)
(define-key ctl-x-map "C" 'ws-copyandset-block)
(define-key ctl-x-map "\C-v" 'ws-copy-block)   ;\C-c
(define-key ctl-x-map "\C-d" 'save-buffers-kill-emacs)
(define-key ctl-x-map "f" 'switch-to-buffer)
(define-key ctl-x-map "\C-f" 'find-file)
(define-key ctl-x-map "\C-l" 'transpose-chars) ;\C-g
(define-key ctl-x-map "\C-h" 'ws-show-markers)
(define-key ctl-x-map "\C-i" 'ws-indent-block)
(define-key ctl-x-map "\C-j" 'help-for-help)
(define-key ctl-x-map "\C-x" 'ws-set-end-of-block)
(define-key ctl-x-map "\C-g" 'ws-downcase-block) ;\C-l
(define-key ctl-x-map "\C-m" 'open-line)
(define-key ctl-x-map "\C-n" 'delete-blank-lines)
(define-key ctl-x-map "\C-o" 'ws-indent-rigidly)
(define-key ctl-x-map "\C-p" 'ws-print-block)
(define-key ctl-x-map "\C-q" 'kill-emacs)
(define-key ctl-x-map "r" 'insert-buffer)
(define-key ctl-x-map "\C-r" 'insert-file)
(define-key ctl-x-map "s" 'save-some-buffers)
(define-key ctl-x-map "\C-s" 'save-buffer)
(define-key ctl-x-map "\C-t" 'ws-mark-word)
(define-key ctl-x-map "\C-u" 'ws-upcase-block)
(define-key ctl-x-map "\C-c" 'ws-move-block)   ;\C-v
(define-key ctl-x-map "\C-w" 'ws-write-block)
(define-key ctl-x-map "\C-k" 'mule-prefix)
(define-key ctl-x-map "\C-y" 'ws-kill-block)
(define-key ctl-x-map "\C-z" 'suspend-emacs)

(define-key ctl-x-map " "    'set-mark-command)
(define-key ctl-x-map "b"    'buffer-menu)
(define-key ctl-x-map "d"    'dired)
(define-key ctl-x-map "k"    'kill-buffer)
(define-key ctl-x-map "V"    'view-file)
(define-key ctl-x-map "w"    'write-file)
(define-key ctl-x-map "_"    'shrink-window)
(define-key ctl-x-map "!"    'shell)
(define-key ctl-x-map "c"    'ws-copy-rectangle)
(define-key ctl-x-map "V"    'ws-copyandset-rectangle)
(define-key ctl-x-map "i"    'ws-open-rectangle)
(define-key ctl-x-map "g"    'ws-clear-rectangle)
(define-key ctl-x-map "v"    'ws-move-rectangle)
(define-key ctl-x-map "y"    'ws-kill-rectangle)
(define-key ctl-x-map ";"    'ws-comment-region)
(define-key ctl-x-map ":"    'ws-uncomment-region)
(define-key ctl-x-map "t"    'toggle-truncate-lines)

  ;; wordstar-C-q-map

(setq ctrl-Q-map (make-keymap))
(fset 'ctrl-Q-prefix ctrl-Q-map)
(defvar ctrl-Q-map nil
  "Ctrl-Q-map maps the key command following ctrl-Q.")

(define-key ctrl-Q-map " " ())
(define-key ctrl-Q-map "0" 'ws-find-marker-0)
(define-key ctrl-Q-map "1" 'ws-find-marker-1)
(define-key ctrl-Q-map "2" 'ws-find-marker-2)
(define-key ctrl-Q-map "3" 'ws-find-marker-3)
(define-key ctrl-Q-map "4" 'ws-find-marker-4)
(define-key ctrl-Q-map "5" 'ws-find-marker-5)
(define-key ctrl-Q-map "6" 'ws-find-marker-6)
(define-key ctrl-Q-map "7" 'ws-find-marker-7)
(define-key ctrl-Q-map "8" 'ws-find-marker-8)
(define-key ctrl-Q-map "9" 'ws-find-marker-9)
(define-key ctrl-Q-map "a" 'ws-query-replace)
(define-key ctrl-Q-map "\C-a" 'isearch-backward)
(define-key ctrl-Q-map "\C-b" 'ws-goto-begining-of-block)
(define-key ctrl-Q-map "\C-v" 'end-of-buffer)                 ;\C-c
(define-key ctrl-Q-map "\C-d" 'end-of-line)
(define-key ctrl-Q-map "\C-e" 'ws-backward-paragraph)
(define-key ctrl-Q-map "f" 'ws-search)
(define-key ctrl-Q-map "\C-f" 'isearch-forward)
(define-key ctrl-Q-map "\C-l" 'transpose-lines)     ;\C-g
(define-key ctrl-Q-map "\C-h" 'ws-kill-bol)
(define-key ctrl-Q-map "\C-j" 'goto-line)
(define-key ctrl-Q-map "\C-x" 'ws-goto-end-of-block)
(define-key ctrl-Q-map "l" 'ws-repeat-search)
(define-key ctrl-Q-map "\C-g" 'recenter) ;\C-l
(define-key ctrl-Q-map "\C-n" 'next-error)
(define-key ctrl-Q-map "\C-o" 'occur)
(define-key ctrl-Q-map "p" 'ws-last-cursorp)
(define-key ctrl-Q-map "\C-p" 'exchange-point-and-mark)
(define-key ctrl-Q-map "\C-q" 'query-replace)
(define-key ctrl-Q-map "\C-r" 'beginning-of-buffer)
(define-key ctrl-Q-map "\C-s" 'beginning-of-line)
(define-key ctrl-Q-map "\C-t" 'backward-kill-word)
(define-key ctrl-Q-map "\C-u" 'replace-string)
(define-key ctrl-Q-map "\C-c" 'overwrite-mode)                ;\C-v
(define-key ctrl-Q-map "w" 'ws-last-error)
(define-key ctrl-Q-map "\C-w" 'ws-backward-sentence)
(define-key ctrl-Q-map "\C-k" 'ws-forward-paragraph)
(define-key ctrl-Q-map "\C-y" 'ws-kill-eol)
(define-key ctrl-Q-map "\C-z" 'ws-forward-sentence)

(define-key ctrl-Q-map "b"    'list-buffers)
(define-key ctrl-Q-map "d"    'list-directory)
(define-key ctrl-Q-map "j"    'what-line)
(define-key ctrl-Q-map "r"    'toggle-read-only)
(define-key ctrl-Q-map "f"    'set-fill-column)

  ;; C-o-map
  
(setq ctrl-O-map (make-keymap))
(fset 'ctrl-O-prefix ctrl-O-map)
(defvar ctrl-O-map nil
  "Ctrl-O-map maps the key command following ctrl-o.")

(define-key ctrl-O-map "0" 'delete-window)
(define-key ctrl-O-map "1" 'delete-other-windows)
(define-key ctrl-O-map "2" 'split-window-vertically)
(define-key ctrl-O-map "3" 'split-window-horizontally)
(define-key ctrl-O-map "4" 'ctl-x-4-prefix)
(define-key ctrl-O-map "5" 'ctl-x-5-prefix)
(define-key ctrl-O-map "\C-a" 'isearch-forward-regexp)
(define-key ctrl-O-map "\C-b" 'auto-fill-mode)
(define-key ctrl-O-map "\C-v" 'scroll-other-window)           ;\C-c
(define-key ctrl-O-map "\C-d" 'scroll-right)
(define-key ctrl-O-map "\C-f" 'isearch-backward-regexp)
(define-key ctrl-O-map "\C-i" 'ws-center-line)
(define-key ctrl-O-map "\C-j" 'justify-current-line)
(define-key ctrl-O-map "o" 'other-window)
(define-key ctrl-O-map "\C-o" 'other-window)
(define-key ctrl-O-map "\C-q" 'query-replace-regexp)
(define-key ctrl-O-map "\C-r" 'scroll-down-other-window)
(define-key ctrl-O-map "\C-s" 'scroll-left)
(define-key ctrl-O-map "\C-u" 'replace-regexp)
(define-key ctrl-O-map "\C-c" 'scroll-up)                     ;\C-v
(define-key ctrl-O-map "\C-w" 'scroll-one-line-down-other-window)
(define-key ctrl-O-map "\C-z" 'scroll-one-line-up-other-window)
(define-key ctrl-O-map "a"    'add-change-log-entry-other-window)
(define-key ctrl-O-map "b"    'switch-to-buffer-other-window)
(define-key ctrl-O-map "d"    'dired-other-window)
(define-key ctrl-O-map "f"    'find-file-other-window)
(define-key ctrl-O-map "m"    'mail-other-window)
(define-key ctrl-O-map "."    'find-tag-other-window)

(setq shell-mode-hook 
      '(lambda ()
	 (define-key shell-mode-map "\C-c\C-v" 'interrupt-shell-subjob) ;\C-v\C-c
	 ))

(setq buffer-menu-mode-hook 
      '(lambda ()
	 (define-key Buffer-menu-mode-map "\C-d" nil)
	 (define-key Buffer-menu-mode-map "\C-x" nil)
	 (define-key Buffer-menu-mode-map "\C-@" 'Buffer-menu-delete-backwards)
	 (define-key Buffer-menu-mode-map "\C-y" 'Buffer-menu-delete)
	 ))

(setq picture-mode-hook 
      '(lambda ()
      (define-key picture-mode-map "\C-f" nil)
      (define-key picture-mode-map "\C-b" nil)
      (define-key picture-mode-map "\C-d" nil)
      (define-key picture-mode-map "\C-c\C-d" nil)
      (define-key picture-mode-map "\C-x" nil)
      (define-key picture-mode-map "\C-o" nil)
      (define-key picture-mode-map "\C-n" nil)
      (define-key picture-mode-map "\C-p" nil)
      (define-key picture-mode-map "\C-e" nil)

      (define-key picture-mode-map "\C-d" 'picture-forward-column)
      (define-key picture-mode-map "\C-s" 'picture-backward-column)
      (define-key picture-mode-map "\C-l" 'picture-clear-column) ;\C-g
      (define-key picture-mode-map "\C-c\C-l" 'delete-char) ;\C-v\C-g
      (define-key picture-mode-map "\C-y" 'picture-clear-line)
      (define-key picture-mode-map "\C-n" 'picture-open-line)
      (define-key picture-mode-map "\C-k" 'picture-move-down)
      (define-key picture-mode-map "\C-e" 'picture-move-up)
      (define-key ctrl-Q-map "\C-d" 'picture-end-of-line)
	 ))

;;(setq fortran-mode-hook 
;;      '(lambda ()
;;	 ))
  
;;(setq f90-mode-hook 
;;      '(lambda ()
;;	 ))
  
(setq c-mode-hook 
      '(lambda ()
	 (define-key c-mode-base-map "\C-d" 'forward-char)
	 (define-key c-mode-base-map "\C-?" 'backward-delete-char-untabify)
	 (c-set-stype "bsd")
	 ))
  
(setq text-mode-hook 
      '(lambda ()
	 (auto-fill-mode 1)
	 ))
	 
(setq TeX-mode-hook 
      '(lambda ()
	 (define-key TeX-mode-map "\C-c\C-c" 'TeX-validate-TeX-buffer) ;\C-v\C-v
	 (setq TeX-show-queue-command "lpstat")
	 (setq TeX-dvi-print-command "dvipr")
	 ))
	 
(setq outline-mode-hook 
      '(lambda ()
	 (define-key outline-mode-map "\C-c\C-h" 'hide-subtree) ;\C-v\C-h
	 ))

(setq dired-mode-hook 
      '(lambda ()
	 (define-key dired-mode-map "\C-d" nil)
	 (define-key dired-mode-map "\C-n" nil)
	 (define-key dired-mode-map "\C-p" nil)
	 (define-key dired-mode-map "\C-l" 'dired-flag-file-deleted) ;\C-g
	 (define-key dired-mode-map "\C-k" 'dired-next-line)
	 (define-key dired-mode-map "\C-e" 'dired-previous-line)
	 ))

(setq view-hook 
      '(lambda ()
	 (define-key view-mode-map "\C-c" 'View-undefined)
	 (define-key view-mode-map "\C-u" 'View-undefined)
	 (define-key view-mode-map "\ev" 'View-undefined)
	 (define-key view-mode-map "\C-p" 'View-undefined)
	 (define-key view-mode-map "\C-s" 'View-undefined)
	 (define-key view-mode-map "\C-q" nil)

	 (define-key view-mode-map "\C-l" 'exit-recursive-edit) ;\C-g
	 (define-key view-mode-map "\C-n" 'universal-argument)
	 (define-key view-mode-map "\C-x" 'Control-X-prefix)
	 (define-key view-mode-map "\C-r" 'View-scroll-lines-backward)
	 (define-key view-mode-map "\C-v" 'View-scroll-lines-forward) ;\C-c
	 (define-key view-mode-map "\C-q\C-g" 'recenter) ; \C-q\C-l
	 (define-key view-mode-map "\C-e" 'View-previous-line)
	 (define-key view-mode-map "\C-k" 'View-next-line)
	 (define-key view-mode-map "\C-q\C-f" 'isearch-forward)
	 (define-key view-mode-map "\C-q\C-a" 'isearch-backward)
	 ))	 

(defun scroll-down-other-window (&optional arg)
  "Scroll down other window."
  (interactive "P")
  (if arg 
      (scroll-other-window (* -1 arg))
    (scroll-other-window (* -2 (window-height)))))

(defun scroll-one-line-up-other-window (&optional arg)
  "Scroll up other window (forward in the text) one line (or N lines)."
  (interactive "p")
  (scroll-other-window(or arg 1)))

(defun scroll-one-line-down-other-window (&optional arg)
  "Scroll down other window (backward in the text) one line (or N)."
  (interactive "p")
  (scroll-down-other-window(or arg 1)))

(defun shrink-window (arg)
  "Shrink window."
  (interactive "p")
  (enlarge-window (* -1 arg)))

(defun ws-forward-sentence ()
  (interactive)
  (if (eq major-mode 'fortran-mode)
      (fortran-next-statement)
    (if (eq major-mode 'f90-mode)
	(f90-next-statement)
      (forward-sentence))))

(defun ws-backward-sentence ()
  (interactive)
  (if (eq major-mode 'f90-mode)
      (f90-previous-statement)
    (if (eq major-mode 'fortran-mode)
	(fortran-previous-statement)
      (backward-sentence))))

(defun ws-forward-paragraph ()
  (interactive)
  (if (eq major-mode 'fortran-mode)
      (fortran-end-of-subprogram)
    (if (eq major-mode 'f90-mode)
	(f90-end-of-subprogram)
      (if (eq major-mode 'c-mode)
	  (end-of-defun)
	(forward-paragraph)))))

(defun ws-backward-paragraph ()
  (interactive)
  (if (eq major-mode 'fortran-mode)
      (fortran-beginning-of-subprogram)
    (if (eq major-mode 'f90-mode)
	(f90-beginning-of-subprogram)
      (if (eq major-mode 'c-mode)
	  (beginning-of-defun)
	(backward-paragraph)))))

(defun ws-mark-paragraph ()
  (interactive)
  (if (eq major-mode 'fortran-mode)
      (fortran-mark-subprogram)
    (if (eq major-mode 'f90-mode)
	(f90-mark-subprogram)
      (if (eq major-mode 'c-mode)
	  (c-mark-function)
	(mark-paragraph)))))

(defun ws-comment-region ()
  (interactive)
  (ws-set-block)
  (if (eq major-mode 'fortran-mode)
      (fortran-comment-region (mark) (point) nil)
    (if (eq major-mode 'f90-mode)
        (f90-comment-region (mark) (point) nil)))
  (ws-off-block))

(defun ws-uncomment-region ()
  (interactive)
  (if (eq major-mode 'fortran-mode)
      (fortran-comment-region (mark) (point) t)
    (if (eq major-mode 'f90-mode)
        (f90-comment-region (mark) (point) t)))
  (ws-set-block))

(make-empty-face 'color-mode-face-a)
(modify-face 'color-mode-face-a "green" nil nil nil nil nil nil)
(point-to-register ?\[ )
(point-to-register ?\] )

(defun ws-set-block ()
  (point-to-register ?\| )
  (register-to-point ?\[ )
  (push-mark (point) nil t)
  (register-to-point ?\] ))

(defun ws-off-block ()
  (pop-mark)
  (register-to-point ?\| ))

(defun ws-put-block-color ()
  (interactive)
  (ws-set-block)
;;  (put-text-property (mark) (point) 'face 'color-mode-face-a)
;;  (put-text-property (mark) (point) 'front-sticky t)
;;  (put-text-property (max (point-min) (1-(point))) (point) 'rear-nonsticky t)
  (ws-off-block))

(defun ws-remove-block-color ()
  (interactive)
  (ws-set-block)
;;  (remove-text-properties (mark) (point) '(face nil))
;;  (remove-text-properties (mark) (min (point-max) (1+(mark))) '(front-sticky nil))
;;  (remove-text-properties (max (point-min) (1-(point))) (point) '(rear-nonsticky nil))
  (ws-off-block))


(defun ws-set-begining-of-block ()
  (interactive)
  (ws-remove-block-color)
  (point-to-register ?\[ )
  (ws-put-block-color)
  (message "Block begin marker set"))

(defun ws-set-end-of-block ()
  (interactive)
  (ws-remove-block-color)
  (point-to-register ?\] )
  (ws-put-block-color)
  (message "Block end marker set"))

(defun ws-goto-begining-of-block ()
  (interactive)
  (register-to-point ?\[ ))

(defun ws-goto-end-of-block ()
  (interactive)
  (register-to-point ?\] ))

(defun ws-write-block ()
  (interactive)
  (ws-set-block)
  (write-region (mark) (point) (read-file-name "Write block to file: "))
  (ws-off-block))

(defun ws-kill-block ()
  (interactive "*")
  (ws-remove-block-color)
  (ws-set-block)
  (kill-region (mark) (point))
  (point-to-register ?\[ )
  (point-to-register ?\] )
  (ws-off-block)
  (ws-put-block-color))

(defun ws-copy-block ()
  (interactive "*")
  (ws-remove-block-color)
  (ws-set-block)
  (setq last-command nil)
  (copy-region-as-kill (mark) (point))
  (ws-off-block)
  (yank)
  (ws-put-block-color))

(defun ws-copyandset-block ()
  (interactive "*")
  (ws-remove-block-color)
  (ws-set-block)
  (setq last-command nil)
  (copy-region-as-kill (mark) (point))
  (ws-off-block)
  (point-to-register ?\[ )
  (yank)
  (point-to-register ?\] )
  (ws-put-block-color))

(defun ws-move-block ()
  (interactive "*")
  (ws-remove-block-color)
  (ws-set-block)
  (setq last-command nil)
  (kill-region (mark) (point))
  (ws-off-block)
  (point-to-register ?\[ )
  (yank)
  (point-to-register ?\] )
  (ws-put-block-color))

(defun ws-show-markers ()
  (interactive)
  (point-to-register ?\| )
  (register-to-point ?\[ )
  (message "Block begin marker")
  (sit-for 1)
  (register-to-point ?\] )
  (message "Block end marker")
  (sit-for 1)
  (register-to-point ?\| )
  (message ""))

(defun ws-indent-block ()
  (interactive "*")
  (ws-set-block)
  (indent-region (mark) (point) nil)
  (ws-off-block))

(defun ws-indent-rigidly ()
  (interactive "*")
  (ws-set-block)
  (indent-rigidly (mark) (point) (string-to-int (read-from-minibuffer "How many columns: ")))
  (ws-off-block))

(defun ws-downcase-block ()
  (interactive "*")
  (ws-set-block)
  (downcase-region (mark) (point))
  (ws-off-block))

(defun ws-upcase-block ()
  (interactive "*")
  (ws-set-block)
  (upcase-region (mark) (point))
  (ws-off-block))

(defun ws-append-to-file ()
  (interactive)
  (ws-set-block)
  (append-to-file (mark) (point) (read-file-name "Append block to file: "))
  (ws-off-block))

(defun ws-append-to-buffer ()
  (interactive)
  (ws-set-block)
  (append-to-buffer (mark) (point) (read-buffer "Append block to buffer: "))
  (ws-off-block))

(defun ws-move-rectangle ()
  (interactive "*")
  (ws-remove-block-color)
  (ws-set-block)
  (kill-rectangle (mark) (point))
  (ws-off-block)
  (point-to-register ?\[ )
  (yank-rectangle)
  (point-to-register ?\] )
  (ws-put-block-color))

(defun ws-copy-rectangle ()
  (interactive "*")
  (ws-remove-block-color)
  (ws-set-block)
  (setq killed-rectangle (extract-rectangle (mark) (point)))
  (register-to-point ?\| )
  (ws-off-block)
  (yank-rectangle)
  (ws-put-block-color))

(defun ws-copyandset-rectangle ()
  (interactive "*")
  (ws-remove-block-color)
  (ws-set-block)
  (setq killed-rectangle (extract-rectangle (mark) (point)))
  (register-to-point ?\| )
  (ws-off-block)
  (point-to-register ?\[ )
  (yank-rectangle)
  (point-to-register ?\] )
  (ws-put-block-color))

(defun ws-open-rectangle ()
  (interactive "*")
  (ws-remove-block-color)
  (ws-set-block)
  (open-rectangle (mark) (point))
  (ws-off-block)
  (ws-put-block-color))

(defun ws-clear-rectangle ()
  (interactive "*")
  (ws-remove-block-color)
  (ws-set-block)
  (clear-rectangle (mark) (point))
  (ws-off-block)
  (ws-put-block-color))

(defun ws-kill-rectangle ()
  (interactive "*")
  (ws-remove-block-color)
  (ws-set-block)
  (kill-rectangle (mark) (point))
  (point-to-register ?\[ )
  (point-to-register ?\] )
  (ws-off-block)
  (ws-put-block-color))

(defun scroll-down-line ()
  (interactive)
  (scroll-down 1))

(defun scroll-up-line ()
  (interactive)
  (scroll-up 1))

;;; ws-mode.el --- WordStar emulation mode for GNU Emacs

;; Copyright (C) 1991 Free Software Foundation, Inc.

;; Author: Juergen Nickelsen <nickel@cs.tu-berlin.de>
;; Version: 0.7
;; Keywords: emulations

;; This file is part of GNU Emacs.

;; GNU Emacs is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.

;; GNU Emacs is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with GNU Emacs; see the file COPYING.  If not, write to
;; the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.

;;; Commentary:

;; This emulates WordStar, with a major mode.

;;;;;;;;;;;
;; wordstar special variables:

(defvar ws-marker-0 nil "Position marker 0 in WordStar mode.")
(defvar ws-marker-1 nil "Position marker 1 in WordStar mode.")
(defvar ws-marker-2 nil "Position marker 2 in WordStar mode.")
(defvar ws-marker-3 nil "Position marker 3 in WordStar mode.")
(defvar ws-marker-4 nil "Position marker 4 in WordStar mode.")
(defvar ws-marker-5 nil "Position marker 5 in WordStar mode.")
(defvar ws-marker-6 nil "Position marker 6 in WordStar mode.")
(defvar ws-marker-7 nil "Position marker 7 in WordStar mode.")
(defvar ws-marker-8 nil "Position marker 8 in WordStar mode.")
(defvar ws-marker-9 nil "Position marker 9 in WordStar mode.")

(defvar ws-search-string nil "String of last search in WordStar mode.")
(defvar ws-search-direction t 
  "Direction of last search in WordStar mode. T if forward, NIL if backward.")

(defvar ws-last-cursorposition nil 
  "Position before last search etc. in WordStar mode.")

(defvar ws-last-errormessage nil 
  "Last error message issued by a WordStar mode function.")

;;;;;;;;;;;
;; wordstar special functions:

(defun ws-center-line ()
  "Center the line point is on, within the width specified by `fill-column'.
This means adjusting the indentation to match
the distance between the end of the text and `fill-column'."
  (interactive)
  (save-excursion
    (let (line-length)
      (beginning-of-line)
      (delete-horizontal-space)
      (end-of-line)
      (delete-horizontal-space)
      (setq line-length (current-column))
      (beginning-of-line)
      (indent-to 
       (+ left-margin 
	  (/ (- fill-column left-margin line-length) 2))))))

(defun ws-error (string)
  "Report error of a WordStar special function. Error message is saved
in ws-last-errormessage for recovery with C-q w."
  (setq ws-last-errormessage string)
  (error string))

(defun ws-set-marker-0 ()
  "In WordStar mode: Set marker 0 to current cursor position."
  (interactive)
  (setq ws-marker-0 (point-marker))
  (message "Marker 0 set"))

(defun ws-set-marker-1 ()
  "In WordStar mode: Set marker 1 to current cursor position."
  (interactive)
  (setq ws-marker-1 (point-marker))
  (message "Marker 1 set"))

(defun ws-set-marker-2 ()
  "In WordStar mode: Set marker 2 to current cursor position."
  (interactive)
  (setq ws-marker-2 (point-marker))
  (message "Marker 2 set"))

(defun ws-set-marker-3 ()
  "In WordStar mode: Set marker 3 to current cursor position."
  (interactive)
  (setq ws-marker-3 (point-marker))
  (message "Marker 3 set"))

(defun ws-set-marker-4 ()
  "In WordStar mode: Set marker 4 to current cursor position."
  (interactive)
  (setq ws-marker-4 (point-marker))
  (message "Marker 4 set"))

(defun ws-set-marker-5 ()
  "In WordStar mode: Set marker 5 to current cursor position."
  (interactive)
  (setq ws-marker-5 (point-marker))
  (message "Marker 5 set"))

(defun ws-set-marker-6 ()
  "In WordStar mode: Set marker 6 to current cursor position."
  (interactive)
  (setq ws-marker-6 (point-marker))
  (message "Marker 6 set"))

(defun ws-set-marker-7 ()
  "In WordStar mode: Set marker 7 to current cursor position."
  (interactive)
  (setq ws-marker-7 (point-marker))
  (message "Marker 7 set"))

(defun ws-set-marker-8 ()
  "In WordStar mode: Set marker 8 to current cursor position."
  (interactive)
  (setq ws-marker-8 (point-marker))
  (message "Marker 8 set"))

(defun ws-set-marker-9 ()
  "In WordStar mode: Set marker 9 to current cursor position."
  (interactive)
  (setq ws-marker-9 (point-marker))
  (message "Marker 9 set"))

(defun ws-find-marker-0 ()
  "In WordStar mode: Go to marker 0."
  (interactive)
  (if ws-marker-0
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-marker-0))
    (ws-error "Marker 0 not set")))

(defun ws-find-marker-1 ()
  "In WordStar mode: Go to marker 1."
  (interactive)
  (if ws-marker-1
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-marker-1))
    (ws-error "Marker 1 not set")))

(defun ws-find-marker-2 ()
  "In WordStar mode: Go to marker 2."
  (interactive)
  (if ws-marker-2
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-marker-2))
    (ws-error "Marker 2 not set")))

(defun ws-find-marker-3 ()
  "In WordStar mode: Go to marker 3."
  (interactive)
  (if ws-marker-3
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-marker-3))
    (ws-error "Marker 3 not set")))

(defun ws-find-marker-4 ()
  "In WordStar mode: Go to marker 4."
  (interactive)
  (if ws-marker-4
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-marker-4))
    (ws-error "Marker 4 not set")))

(defun ws-find-marker-5 ()
  "In WordStar mode: Go to marker 5."
  (interactive)
  (if ws-marker-5
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-marker-5))
    (ws-error "Marker 5 not set")))

(defun ws-find-marker-6 ()
  "In WordStar mode: Go to marker 6."
  (interactive)
  (if ws-marker-6
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-marker-6))
    (ws-error "Marker 6 not set")))

(defun ws-find-marker-7 ()
  "In WordStar mode: Go to marker 7."
  (interactive)
  (if ws-marker-7
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-marker-7))
    (ws-error "Marker 7 not set")))

(defun ws-find-marker-8 ()
  "In WordStar mode: Go to marker 8."
  (interactive)
  (if ws-marker-8
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-marker-8))
    (ws-error "Marker 8 not set")))

(defun ws-find-marker-9 ()
  "In WordStar mode: Go to marker 9."
  (interactive)
  (if ws-marker-9
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-marker-9))
    (ws-error "Marker 9 not set")))

(defun ws-mark-word ()
  "In WordStar mode: Mark current word as block."
  (interactive)
  (save-excursion
    (forward-word 1)
    (sit-for 1)
    (ws-set-end-of-block)
    (forward-word -1)
    (sit-for 1)
    (ws-set-begining-of-block)))

(defun ws-goto-last-cursorposition ()
  "In WordStar mode: "
  (interactive)
  (if ws-last-cursorposition
      (let ()
	(setq ws-last-cursorposition (point-marker))
	(goto-char ws-last-cursorposition))
    (ws-error "No last cursor position available.")))

(defun ws-last-error ()
  "In WordStar mode: repeat last error message.
This will only work for errors raised by WordStar mode functions."
  (interactive)
  (if ws-last-errormessage
      (message ws-last-errormessage)
    (message "No WordStar error yet.")))

(defun ws-kill-eol ()
  "In WordStar mode: Kill to end of line (like WordStar, not like Emacs)."
  (interactive)
  (let ((p (point)))
    (end-of-line)
    (kill-region p (point))))

(defun ws-kill-bol ()
  "In WordStar mode: Kill to beginning of line 
\(like WordStar, not like Emacs)."
  (interactive)
  (let ((p (point)))
    (beginning-of-line)
    (kill-region (point) p)))

(defun kill-complete-line ()
  "Kill the complete line."
  (interactive)
  (beginning-of-line)
  (if (eobp) (error "End of buffer"))
  (let ((beg (point)))
    (forward-line 1)
    (kill-region beg (point))))

(defun ws-search (string)
  "In WordStar mode: Search string, remember string for repetition."
  (interactive "sSearch for: ")
  (message "Forward (f) or backward (b)")
  (let ((direction
	 (read-char)))
    (cond ((equal (upcase direction) ?F)
	   (setq ws-search-string string)
	   (setq ws-search-direction t)
	   (setq ws-last-cursorposition (point-marker))
	   (search-forward string))
	  ((equal (upcase direction) ?B)
	   (setq ws-search-string string)
	   (setq ws-search-direction nil)
	   (setq ws-last-cursorposition (point-marker))
	   (search-backward string))
	  (t (keyboard-quit)))))

(defun ws-repeat-search ()
  "In WordStar mode: Repeat last search."
  (interactive)
  (setq ws-last-cursorposition (point-marker))
  (if ws-search-string
      (if ws-search-direction
	  (search-forward ws-search-string)
	(search-backward ws-search-string))
    (ws-error "No search to repeat")))

(defun ws-query-replace (from to)
  "In WordStar mode: Search string, remember string for repetition."
  (interactive "sReplace: 
sWith: " )
  (setq ws-search-string from)
  (setq ws-search-direction t)
  (setq ws-last-cursorposition (point-marker))
  (query-replace from to))

;;; ws-mode.el ends here
