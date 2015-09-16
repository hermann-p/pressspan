(ns pressspan.graph
  (:require [clj-orient.core :as oc]
            [clj-orient.graph :as og]
            [clj-orient.query :as oq]
            [clojure.data.int-map :as im]))

(use '[clojure.test]
     '[taoensso.timbre.profiling])


(declare select-fragment
         add-chromosome
         register-fragment)


(def db-name "memory:test")
(def db-user "admin")
(def db-pass "admin")


(def rid #(:#rid %)) ; workaround for linebreak-bug of :#rid


(defn init-database!
  "Opens a database from global parameters
  db-name, db-user, db-pass.
  Initialises global parameters
  circulars, multis
  Creates required datatypes to store a genome."
  []
  ; Side-effects!
  (def circulars (ref #{}))
  (def multis (ref #{}))
  
  (og/create-vertex-type! :genome)
    (og/create-vertex-type! :chromosome)
    (og/create-vertex-type! :fragment)
    (og/create-edge-type! :splice)
    (og/create-edge-type! :has-chr)
    (og/create-edge-type! :has-frag)
    (oc/create-prop! :chromosome :p5frags :link-map)
    (oc/create-prop! :chromosome :p3frags :link-map)
    (oc/create-prop! :genome :chrs :link-map)
    (oc/save-schema!)
    
    (try
      (oc/with-tx
        (let [root (oc/save! (og/vertex :genome))]
          (og/add-root! :genome root)))
      
      (catch Exception e
        (println e))))


(defnp process-header
  [lines head? process]
  (let [inserter
        (fn [chrs chr-data]
          (let [chr-id (:id chr-data)
                chr-data (assoc chr-data
                           :p3frags (im/int-map)
                           :p5frags (im/int-map)
                           :mp3frags (im/int-map)
                           :mp5frags (im/int-map))]
            (assoc! chrs chr-id chr-data)))

        all-chromosomes
        (reduce inserter (transient {})
                (filter (complement nil?)
                        (for [line lines :while (head? line)]
                          (process line))))]
    (persistent! all-chromosomes)))


(defnp process-data
  [chrs lines get-data add-all head?]
  (let [inserter
        (fn [chrs frag]
          (let [[l5 l3] (if (= :plus (:dir frag))
                          [:p5frags :p3frags]
                          [:mp5frags :mp3frags])
                frag (add-all chrs frag)
                orid (rid frag)
                chr (get chrs (:chr frag))]
            (assoc! chrs (:chr frag) 
                    (assoc chr
                      l5 (into (l5 chr) {(:p5 frag) orid})
                      l3 (into (l3 chr) {(:p3 frag) orid})))))
        all-chrs
        (reduce inserter (transient chrs)
                (for [line (drop-while head? lines)]
                  (get-data line)))]
    (persistent! all-chrs)))


(defn create-genome
  "Connects to an OrientDB database and calls functions to read in a genome file
  Input:
  [file-name] string path to datafile
  [fns]       map of processing functions. Needs
              :head?          predicate to identify header lines
              :header-parser  parse a header line string for a chromosome
              :data-parser    parse a data line for fragment information
              :add-all        composition of insertion- and filter functions
                              or nil to only insert"
  [gname file-name fns & credentials]
  (def db-name gname)
  (println "Processing" file-name)
  (if credentials
    (do
      (println "Using provided database credentials")
      (def db-user (:user credentials))
      (def db-pass (:pass credentials))))
  (println "Requested OrientDB database:" db-name)
  (let [db-type (first (clojure.string/split db-name #":"))]
    (println "OrientDB type:" db-type)
    (if (contains? #{"memory" "local" "plocal"}  db-type)
      (do
        (println "Creating fresh database...")
        (oc/create-db! gname))
      (println "Trying to connect to database...")))
  (println "Opening" db-name)
  (oc/with-db
    (og/open-graph-db! db-name db-user db-pass)
    (init-database!)
    (if-let [file (pressspan.io/lazy-file file-name)]
      (def chromosomes
        (-> (process-header file
                            (:head? fns)
                            (:header-parser fns))
            (process-data file
                          (:data-parser fns)
                          (or (:add-all fns) register-fragment)
                          (:head? fns)))))))


(defnp select-chromosome
  [chr id]
  (if id (get chr id)))


(defnp link-frags
  "link [up-el] => [dn-el]
  [count] increment link counter if non-nil
  returns the link"
  [up-el dn-el & count]
  (if-let [link (first (og/get-links up-el dn-el))]
    (if count (oc/save! (update link :depth inc)))
    (oc/save! (og/link! up-el :splice {:depth 1} dn-el))))


(defnp select-fragment
  "selects and returns a fragment from the database
   [chr] chromosome or chromosome id-string
   [dir] :p3 or :p5 for 3' or 5' end
   [pos] position of 3' or 5' end"
  [chrs chr dir pos strandiness]
  (if-let [chrom (select-chromosome chrs (if (map? chr) (:id chr) chr))] ; select chromosome if neccessary, proceed only if found
    (let [end (if (= strandiness :plus)
                (if (= dir :p5) :p5frags :p3frags)
                (if (= dir :p5) :mp5frags :mp3frags))]
      (if-let [orid (get (end chrom) pos)]
        (try                                                             ; catch exception for REPL tests
          (oc/load orid)
          (catch Exception e
            orid))))))


(defnp register-fragment
  "creates a fragment, registers and links it within the chromosome database
  performs a check before creation to avoid duplicates
  [data] hash-map of fragment properties"
  [chrs {:keys [prev next chr p3 p5 dir]}]    ; dir in {:plus, :minus}, hash-maps {prev, next}
  (oc/with-tx
    (let [frag (or (select-fragment chrs chr :p5 p5 dir)
                   (oc/save! (og/vertex :fragment {:p3 p3, :p5 p5, :chr chr :dir dir})))
          chr (select-chromosome chrs chr)]
      (if next
        (if-let [nfrag (select-fragment chrs (:chr next) :p5 (:p5 next) (:dir next))]
          (link-frags frag nfrag)))
      (if prev
        (if-let [pfrag (select-fragment chrs (:chr prev) :p3 (:p3 prev) (:dir prev))]
          (link-frags pfrag frag :count)))
    (assoc frag :next next :prev prev))))      ; return document and additional info, for composition with filters


(defnp remember-multistrand
  [frag]
  (if-not (= (:chr frag) (:chr (:next frag))) ; multistrand if elements on different strands
    (dosync
     (commute multis conj (:#rid frag)))) frag)


(defnp remember-circular
  [frag]
  (if-let [next-p5 (:p5 (:next frag))]
    (let [is-upstream? (if (= :plus (:dir frag)) < >)] ; define upstream-test
      (if ((every-pred true?)
               (= (:chr frag) (:chr (:next frag)))     ; circular if on same strand and
               (is-upstream? next-p5 (:p5 frag)))      ; next frag starts upstream on same strand
        (dosync
         (commute circulars conj (:#rid frag)))))) frag)


;; Traversing the database requieres casting to and from record-ids
;; as different db fetches of the same element may or may not be unique

(defn- get-neighbours
  "Helper to traverse a graph from the database."
  [el]
  (let [n-coll (for [edge (og/get-edges (oc/load el) :both :splice)]
                 ((juxt #(og/get-vertex % :in)
                        #(og/get-vertex % :out))
                  edge))]
    (into #{} (map :#rid (apply concat n-coll)))))
(defn get-subgraph
  "Get a sequence of frags linked to [start] using a BFS approach
  Needs access to an open database."
  [start]
  (let [walk
        (fn walk [known this]
          (if-let [next-el (first (clojure.set/difference (get-neighbours this) known))]
            (lazy-seq
             (cons this
                   (walk (conj known this) next-el)))))] (walk #{} start)))


(deftest sam-test                ; run tests on manually created mapping file
  (if (oc/db-exists? "memory:test")
    (do (println "\nDeleting database...")
        (oc/delete-db! (og/open-graph-db! "memory:test" "admin" "admin"))))
  (create-genome "memory:test" "test/data/5_out.sam"
                 {:head?         pressspan.saminput/header-line?
                  :header-parser pressspan.saminput/parse-header-line
                  :data-parser   pressspan.saminput/make-frag
                  :add-all (comp remember-circular remember-multistrand register-fragment)}) ; test realignment, multistrand and circular detection
  (oc/with-db (og/open-graph-db! "memory:test" "admin" "admin")
    (is (= 14 (count (oq/native-query :fragment {})))) ; 14 unique fragments, 6 duplicates
    (is (< 0 (og/count-edges)))                        ; links should appear
    (is (= 3 (count (get-subgraph (first @multis)))))  ; first multisplit consists of 3 elements
    (println (clojure.string/join \newline (get-subgraph (first @multis))))) ; look at the links
  (is (= 12  (count @multis))))                        ; there are 12 multistrand junctions

; (run-tests)

(if false ; code to print all created fragments ordered by chromosome
  (oc/with-db (og/open-graph-db! "memory:test" "admin" "admin")
    (for [frag (sort #(compare (:chr %1) (:chr %2)) (oq/native-query :fragment {}))]
      (println (format "%s:%d-%d" (:chr frag) (:p5 frag) (:p3 frag))))))
