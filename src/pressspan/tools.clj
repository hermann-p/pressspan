(ns pressspan.tools)

(defn rev-comp
  "Calculates the reverse complement of a sequence string"
  [s]
  (clojure.string/join "" (for [c (clojure.string/upper-case (reverse s))]
         (case c
           \A \U
           \T \A
           \G \C
           \C \G
           ""))))
