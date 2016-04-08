(ns pressspan.statistics
  (:require [clojure.java.io :refer [writer]]))

(defn add-empty-stats
  "Create zero-initialized statistics row with [len] buckets"
  [input n]
  (let [id (first input)
        len (:len (second input))
        n (if (seq? n) (first n) n)]
    (println "len" (type len) len "n" (type n) n)
    {id
     {:id id
      :bsize (max 1 (int (/ len n)))
      :vals (vec (take n (cycle [0])))}}))

(defn empty-stats
  "Create zero-initialized statistics table for all chromosomes"
  [genome n-buckets]
  (reduce into {} (map #(add-empty-stats % n-buckets)
                       (dissoc genome :frags :links :multis :circulars :custom :nl :nf))))

(defn get-line
  "Retrieve the nth sorted line of a statistics table"
  [n s]
  (for [k (sort (keys s))]
    (get-in s [k :vals n])))

(defn get-bucket-sizes
  "Retrieve a sorted line of bucket sizes for statistics table"
  [s]
  (for [k (sort (keys s))]
    (get-in s [k :bsize])))
    
(defn bucket-size
  "Retrieve bucket size of nth chromosome"
  [root n]
  {:pre [(identity (get root n))
         (get-in root [n :bsize])]}
  (get-in root [n :bsize]))

(defn increase
  "Increase number of counts for element"
  [root stat chr pos dir]
  (dosync
   (let [bs (get-in @stat [chr :bsize])
         pos (int (if (= :plus dir)
                    (/ pos bs)
                    (/ (- (get-in @root [chr :len]) pos) bs)))]
     (commute stat update-in [chr :vals pos] inc))))

(defn write-stat-file
  "Generate a statistics table for a event type from a genome and
  write it to a file. [seed-key] from {:multis, :circulars} defines
  type of events"
  [genome file-name seed-key]
  (let [stats (sort (if (= :multis seed-key) (genome :mstat) (genome :cstat)))]
    (println "Writing statistics:" file-name)
    (with-open [wrtr (writer file-name)]
      (.write wrtr (clojure.string/join "\t" (keys stats)))
      (.write wrtr "\n")
      (.write wrtr (clojure.string/join "\t" (get-bucket-sizes stats)))
      (.write wrtr "\n")
      (let [stats (map second stats)]
        (loop [rem (map :vals stats)]
          (when (seq (first rem))
            (.write wrtr
                    (clojure.string/join
                     "\t"
                     (map first rem)))
            (.write wrtr "\n")
            (recur (map rest rem)))))))
  genome)
