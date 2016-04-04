(ns pressspan.io
  (:use [clojure.java.io :only [reader]]
        clojure.test))

(defn lazy-file
  "Wrapper for memory-efficient lazy file reading"
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


(defn make-dir
  "Create a directory if it doesn't exist yet"
  [path]
  (if-not (.mkdirs (java.io.File. path))
    (println "WARNING: cannot create directory" path
             "- does it exist already?"))); .mkdirs <-> mkdir -p


(deftest io-tests
  (let [lf (lazy-file "tests/data/5_out.sam")]
    (is (complement (nil? lf)))
    ))

; (run-tests)
