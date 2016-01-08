(ns pressspan.statistics
  (:require [clojure.java.io :refer [writer]]))

(defn add-empty-stats [input n]
  (let [id (first input)
        len (:len (second input))]
    {id
     {:id id
      :bsize (int (/ len n))
      :vals (vec (take n (cycle [0])))}}))

(defn empty-stats [genome n-buckets]
  (reduce into {} (map #(add-empty-stats % n-buckets)
                       (dissoc genome :frags :links :multis :circulars :custom :nl :nf))))


(defn get-linked [genome link-id]
  (let [el-id (get-in genome [:links link-id :down])]
    (get-in genome [:frags el-id])))


(defn is-multi? [a b]
  {:pre [(map? a) (map? b)]}
  ((some-fn true?) (not= (:chr a) (:chr b)) (not= (:dir a) (:dir b))))


(defn is-circular? [a b]
  {pre [(map? a) (map? b)]}
  (every? true? [(= (:chr a) (:chr b)) (= (:dir a) (:dir b)) (< (:p3 b) (:p3 a))]))


(defn get-pairs [genome pred? seeds]
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


(defn get-line [n s]
  (for [k (sort (keys s))]
    (get-in s [k :vals n])))

(defn get-bucket-sizes [s]
  (for [k (sort (keys s))]
    (get-in s [k :bsize])))
    
(defn bucket-size [root n]
  {:pre [(identity (get root n))
         (get-in root [n :bsize])]}
  (get-in root [n :bsize]))


(defn write-stat-file
  [genome file-name seed-key n-buckets]
  (println "Writing statistics:" file-name "with" n-buckets "buckets per chromosome...")
  (let [[pred seeds] (if (= :multis seed-key)
                       [is-multi? (:multis genome)]
                       [is-circular? (:circulars genome)])
        stats (reduce
               (fn [s vals]
;;                {:pre [(seq s) (seq vals)]}
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
          (recur (inc i)))))))
