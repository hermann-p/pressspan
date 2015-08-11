(ns pressspan.io)

(use '[clojure.java.io :only [reader]]
     '[clojure.test])

(defn lazy-file
  [file-name]
  (letfn [(helper [rdr]
            (lazy-seq
             (if-let [line (.readLine rdr)]
               (cons line (helper rdr))
               (do (.close rdr) nil))))]
    (try
      (helper (reader file-name))
      (catch Exception e
        (println "Exception while trying to open file" file-name)
        ))))


(deftest io-tests
  (let [lf (lazy-file "/home/hermann/5_out.sam")]
    (is (complement (nil? lf)))
    ))

; (run-tests)
