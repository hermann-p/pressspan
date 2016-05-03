(ns pressspan.visualise
  (:use [clojure.test]
        [clojure.set :only [difference union]])
  (:require [pressspan.graph :as graph]
            [clojure.java.io :refer [writer]]))


(defn- hex
  "Converts integer n to two-digit hexadecimal number"
  [n]
  {:pre [(integer? n) (>= n 0) (<= n 255)]}
  (if (> n 16)
    (format "%x" n)
    (format "0%x" n)))

(defn- chr-color
  "Automatically generates color for nth chromosome"
  [n]
  {:pre [(integer? n)]}
  (nth (cycle (for [r (map (partial * 85) (range 4))
                    g [255 0 128]
                    b (map #(- 255 (* 85 %)) (range 4))]
                (str "#" (hex r) (hex g) (hex b))))
       (* 3 n) 1))


(defn graph->log
  "Creates a tab-seperated CSV string from a graph"
  [root graph id-string]
  (format "%s\t|%s|\t%s\n"
          id-string 
	  (clojure.string/join ","
                               (for [node graph] 
                                 (str (if (= :minus (:dir node)) "-") (:chr node) ":" (:p5 node) "-" (:p3 node))))
          (clojure.string/join ","
                               (sort (set (map :chr graph))))))

(declare color-codes) ;; forward declaration
(defn graph->dot
  "Creates a .dot format string from a graph"
  [root graph id-string]
  (let [put-el
        (fn [el]
          (cons 
           (str "  " (:id el) " [shape=\"rectangle\",border=0,style=\"" 
                (if (= :plus (:dir el)) "filled" "bold") "\",height=0.25,"
                "label=\"" (:chr el) ":" (:p5 el) "-" (:p3 el) "\","
                "color=\"" (chr-color (get color-codes (:chr el))) "\"];\n")
           (for [link (:up el)]
             (str "  " (:id el) "->" (:up link) " [label=\"" (:depth link) "\"];\n"))))
        walk
        (fn walk [els]
          (if (seq els)
            (lazy-seq
             (conj (walk (rest els)) (clojure.string/join (put-el (first els)))))))]
    (str "digraph " id-string " {" \newline 
         "  rankdir=\"LR\";" \newline
         (clojure.string/join (doall (concat (walk graph))))
         "}")))


(defn drop-filter
  "Filter to sort out graphs in which any edge has less than [min-depth] reads"
  [min-depth]
  {:pre [(integer? min-depth)]}
  (fn [graph]
    (every? true?
            (for [node graph]
              (or (empty? (:up node))
                  (every? #(<= min-depth %) (map :depth (:up node))))))))

(defn range-filter
  "Filter all graphs on chromosome chr, strandiness dir, p5 or p3
   between upper and lower"
  [chr dir lower upper]
  {:pre [(string? chr)
         (#{:plus :minus} dir)
         (integer? upper)
         (integer? lower)
         (< lower upper)]}
  (fn [graph]
    (let [candidates (filter #(= chr (:chr %)) graph)
          candidates (filter #(= dir (:dir %)) candidates)]
      (if (empty? candidates) false
      	  (let [matches 
                                        ; test nodes if the read overlaps the boundaries partially or completely
      	        (map #(some-fn (<= lower (:p5 %) upper)
                               (<= lower (:p3 %) upper)
                               ((every-pred (fn [el] (< (:p5 el) lower))
                                            (fn [el] (> (:p3 el) upper))) %))
                     candidates)]
            (not (every? false? matches)))))))


(defn write-files
  "Writes all subgraphs of a chosen category to .dot files"
  ([root type basedir] (write-files root type basedir 1))
  ([root type basedir min-depth & filters]
   (def color-codes                   ; create a lookup-map for colors from chromosome names
     (let [chromosomes (filter string? (keys root))]
       (zipmap chromosomes (range))))
   (let [meaningful?
         (fn [g] (< 1 (count g) 100))   ; graph has at least two and no more than 100 nodes

         filters (flatten (into filters [meaningful?]))
                                        ;		  		filters (filter identity filters)

         passes?
         (fn [g]
           {:pre (every? identity filters)}
           (every? true? ((apply juxt filters) g)))

         graphs (filter passes? (graph/all-subgraphs root (type root) min-depth))
         typestr (name type)]
     (when-not (zero? (count graphs))
       (println "Writing" (count graphs) typestr "graph files to" (str basedir "/" typestr))
       (pressspan.io/make-dir (str basedir "/" typestr))
       (doall
        (map-indexed
         #(spit
           (str basedir "/" typestr "/" typestr "_" %1 ".dot")
           (graph->dot root %2 (str typestr "_" %1)))
         graphs))
       (print "-> appending to log file... ")
       (with-open [wrtr (writer (str basedir "/pressspan.log") :append true)]
         (doall
          (map-indexed
           #(.write wrtr
                    (graph->log root %2 (str typestr "_" %1)))
           graphs)))
       (println "done.")))
   root))
