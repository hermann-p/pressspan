(ns pressspan.statistics
  (:require [clojure.java.io :refer [writer]]))

(defn add-empty-stats
  "Create zero-initialized statistics row with [len] buckets"
  [input n]
  (let [id (first input)
        len (:len (second input))]
    {id
     {:id id
      :bsize (max 1 (int (/ len n)))
      :vals (vec (take n (cycle [0])))}}))

(defn empty-stats
  "Create zero-initialized statistics table for all chromosomes"
  [genome n-buckets]
  (reduce into {} (map #(add-empty-stats % n-buckets)
                       (dissoc genome :frags :links :multis :circulars :custom :nl :nf))))


(defn get-linked
  "Retrieve all downstream linked fragments of a fragment"
  [genome link-id]
  (let [el-id (get-in genome [:links link-id :down])]
    (get-in genome [:frags el-id])))


(defn is-multi?
  "Are fragments a and b linked by multistrand splicing?"
  [a b]
  {:pre [(map? a) (map? b)]}
  ((some-fn true?) (not= (:chr a) (:chr b)) (not= (:dir a) (:dir b))))


(defn is-circular?
  "Are a and b linked by backsplicing?"
  [a b]
  {pre [(map? a) (map? b)]}
  (let [upstream? (if (= :plus (:dir a)) < >)]
    (every? true? [(= (:chr a) (:chr b))
                   (= (:dir a) (:dir b))
                   (upstream? (:p5 b) (:p5 a))])))


(defn get-pairs
  "Get both ends of all splicing of type pred {is-multi?, is-circular?}"
  [genome pred? seeds]
  (->> (pmap
        (fn [s]
          (let [s (get-in genome [:frags s])]
            (for [l (:down s)]
              (let [link (get-linked genome l)]
                (if (pred? s link)
                  [[(:chr s) (:p5 s)]
                   [(:chr link) (:p3 link)]
                   (get-in genome [:links l :depth])])))))
        seeds)
       (filter #(every? identity [(seq %) (identity (first %))]))))

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


(defn write-stat-file
  "Generate a statistics table for a event type from a genome and
  write it to a file. [seed-key] from {:multis, :circulars} defines
  type of events"
  [genome file-name seed-key n-buckets]
  (println "Writing statistics:" file-name "with" n-buckets "buckets per chromosome...")
  (let [[pred seeds] (if (= :multis seed-key)
                       [is-multi? (:multis genome)]
                       [is-circular? (:circulars genome)])
        stats (reduce
               (fn [s vals]
                (if-not (seq vals)
                 s
                 (let [[[id-a pa] [id-b pb] l] (first vals)
                       end-a (min (int (/ pa (bucket-size s id-a)))
                                  (dec n-buckets))
                       start-b (min (int (/ pb (bucket-size s id-b)))
                                    (dec n-buckets))]
                   (-> (update-in s [id-a :vals end-a] #(+ % l))
                       (update-in [id-b :vals start-b] #(+ % l))))))
               (empty-stats genome n-buckets)
               (get-pairs genome pred seeds))]
    (with-open [wrtr (writer file-name)]
      (.write wrtr (clojure.string/join "\t" (sort (keys stats))))
      (.write wrtr "\n")
      (.write wrtr (clojure.string/join "\t" (get-bucket-sizes stats)))
      (.write wrtr "\n")
      (loop [i 0]
        (when (< i n-buckets)
          (.write wrtr
            (clojure.string/join
             "\t"
             (get-line i stats)))
          (.write wrtr "\n")
          (recur (inc i))))))
    genome)
