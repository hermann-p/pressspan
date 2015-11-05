(ns pressspan.visualise
  (:use [clojure.test]
        [clojure.set :only [difference union]])
  (:require [pressspan.graph :as graph]
            [clojure.java.io :refer [writer]]))


(defn assign-layers
  "returns graph element with assigned layers.
   root  - map create by pressspan.graph/create-genome
   graph - seq of element indices"
  [root graph]
  (letfn
    [(up-known? [up-known verts el-id]
                (let [el (get verts el-id)
                      links (map #(get-in root [:links %]) (:up el))
                      this-up (set (map :up links))]
                  (empty? (difference this-up up-known))))    ; el has no unknown upstream links
     
     (find-next [free verts up-els]
                (first (filter (partial up-known? up-els verts) free)))
     
     (find-here [free up-els els-here vert-map layered-graph layer]
                (if (seq free)
                  (if-let [v (find-next free vert-map up-els)]
                               
                    (find-here                            ; map found el to current layer
                      (disj free v)                            ; remove it from available els
                      up-els                                   ; leave upper layer elements alone
                      (conj els-here v)                        ; store el in els of this layer
                      (assoc-in vert-map [v :layer] layer)     ; assign layer number to el
                      (update-in layered-graph [layer] conj v) ; add element to layered graph
                      layer)                                   ; keep layer number
                    (find-here                            ; step up one layer
                      free                                     ; keep free elements
                      (union up-els els-here)                  ; save elements in "known from upper layer"
                      #{}                                      ; empty set of els in new layer
                      vert-map                                 ; keep vertices
                      layered-graph                            ; keep layered graph
                      (inc layer)))                            ; increase layer number
                  layered-graph))]
    (let [vertices (mapv #(get-in root [:frags %]) graph)
          vert-map (zipmap graph vertices)]
      (trampoline find-here (set graph) #{} #{} vert-map {} 0 ))))


(defn- hex
  "Converts integer n to two-digit hexadecimal number"
  [n]
  {:pre [(integer? n) (>= n 0) (<= n 255)]}
  (if (> n 16)
    (format "%x" n)
    (format "0%x" n)))

(defn- chr-color [n]
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

(declare color-codes)
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
  ; if lower <= x <= upper, then one of [lower-x, upper-x] is <= 0 and one >= 0
  (fn [graph]
    (let [candidates (filter #(= chr (:chr %)) graph)
          candidates (filter #(= dir (:dir %)) candidates)]
        (println "checking range: " graph)
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
  	(println "Applied filters:" \newline (clojure.string/join "\n " filters))
    (println "Writing" (count graphs) typestr "graph files to" (str basedir "/" typestr))
    (pressspan.io/make-dir (str basedir "/" typestr))
    (doall
      (map-indexed
        #(spit
            (str basedir "/" typestr "/" typestr "_" %1 ".dot")
            (graph->dot root %2 (str typestr "_" %1)))
        graphs))
    (with-open [wrtr (writer (str basedir "/pressspan.log") :append true)]
      (doall
        (map-indexed
          #(.write wrtr
            (graph->log root %2 (str typestr "_" %1)))
          graphs)))
  root)))


(deftest graph-test
  (let [filename "test/data/5_out.sam"
        funs {:head? pressspan.saminput/header-line?
              :header-parser pressspan.saminput/parse-header-line
              :data-parser pressspan.saminput/make-frag
              :add-all [pressspan.graph/remember-multistrand pressspan.graph/register-fragment]}
        genome (pressspan.graph/create-genome filename funs)
        multi (first (:multis genome))
        subgraph (pressspan.graph/get-subgraph genome multi)]
    (is (true? ((drop-filter 2) [{:up [{:depth 2} {:depth 5}]} {:up [{:depth 3}]} {:down nil}])))
    (is (false? ((drop-filter 3) [{:up [{:depth 2} {:depth 5}]} {:up [{:depth 3}]} {:down nil}])))
    (write-files genome :multis "test/output" 1)
    (println (graph->dot genome subgraph "testgraph"))
    (println (graph->log genome subgraph "testgraph"))))