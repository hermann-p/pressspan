(ns pressspan.graph
  (:require [clj-orient.core :as oc]
            [clj-orient.graph :as og]
            [clj-orient.query :as oq]
            [clojure.data.int-map :as im]))
;  (:require
;    [taoensso.timbre :as timbre
;      :refer (log  trace  debug  info  warn  error  fatal  report
;              logf tracef debugf infof warnf errorf fatalf reportf
;              spy get-env log-env)]
;    [taoensso.timbre.profiling :as profiling
;      :refer (pspy pspy* profile defnp p p*)])) ; hm... seems to cast int-map to map?

(use '[clojure.test]
     '[taoensso.timbre.profiling])

(declare graph-test-fn
         select-fragment
         add-chromosome
         register-fragment)


(def db-name "memory:test")
(def db-user "admin")
(def db-pass "admin")


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
  (doseq [line lines :while (head? line)]
    (if-let [info (process line)]
      (add-chromosome (:id info) (:len info)))))


(defnp process-data
  [lines get-data add-all head?]
  (doseq [line (drop-while head? lines)]
    (add-all (get-data line))))


(defn create-genome
  "Cennects to an OrientDB database and calls functions to read in a genome file
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
  (if credentials
    (do
      (println "Using provided database credentials")
      (def db-user (:user credentials))
      (def db-pass (:pass credentials))))
  (println "Requested OrientDB database:" db-name)
  (let [db-type (first (clojure.string/split db-name #":"))]
    (println "OrientDB type:" db-type)
    (if ((some-fn true?)
         (= db-type "memory")
         (= db-type "local")
         (= db-type "plocal"))
      (do
        (println "Creating fresh database...")
        (oc/create-db! gname))
      (println "Trying to connect to database...")))
  (println "Opening" db-name)
  (oc/with-db
    (og/open-graph-db! db-name db-user db-pass)
    (init-database!)
    (if-let [file (pressspan.io/lazy-file file-name)]
      (do
        (process-header file (:head? fns) (:header-parser fns))
        (process-data file (:data-parser fns) (or (:add-all fns) register-fragment) (:head? fns))))))


(defnp add-chromosome
  "creates an empty chromosome with lookup-id [id] and a length of [len] basepairs,
   registers and links it with the genome graph"
  [id len]
  (let [chr-data {:id id
                  :len len
                  :p5frags (im/int-map)   ; 5' +strand links
                  :p3frags (im/int-map)   ; 3' +strand links
                  :mp5frags (im/int-map)  ; 5' -strand links
                  :mp3frags (im/int-map)} ; 3' -strand links
        chr (oc/with-tx (oc/save! (og/vertex :chromosome chr-data)))
        root (og/get-root :genome)]
    (oc/with-tx
      (oc/save! (assoc root :chrs (assoc (:chrs root) id (:#rid chr)))))))

    
(defnp select-chromosome
  "selects a chromosome by [id]-string from the database and returns it"
  [id]
  (if id
    (try
      (oc/load (get (:chrs (og/get-root :genome)) id))            ; throws if no such chromosome
      (catch Exception e
        (println "select-chromosome - No such chromosome:" id)
        nil))))


(defnp link-frags
  "link [up-el] => [dn-el]
  [count] increment link counter if non-nil
  returns the link"
  [up-el dn-el & count]
;  (println "Linkin" (dissoc up-el :in :out) "=>" (dissoc dn-el :in :out))
  (oc/with-tx
    (if-let [link (first (og/get-links up-el dn-el))]
      (if count (oc/save! (update link :depth inc)))
      (oc/save! (og/link! up-el :splice {:depth 1} dn-el)))))


;(defnp link-with-chr
;  [frag link pos chr]
;  (let [newlink (assoc chr link (p :newlink (assoc (link chr) pos (:#rid frag))))] (p :savelink (oc/save! newlink))))

(defnp link-with-chr
  [frag link pos chr]
  (oc/with-tx
    (oc/save! (update-in chr [link] into {pos (:#rid frag)}))))


(defnp select-fragment
  "selects and returns a fragment from the database
   [chr] chromosome or chromosome id-string
   [dir] :p3 or :p5 for 3' or 5' end
   [pos] position of 3' or 5' end"
  [chr dir pos strandiness]
  (if-let [chrom (select-chromosome (if (map? chr) (:id chr) chr))] ; select chromosome if neccessary, proceed only if found
    (let [end (if (= strandiness :plus)
                (if (= dir :p5) :p5frags :p3frags)
                (if (= dir :p5) :mp5frags :mp3frags))
          pos (str pos)]
      (if-let [orid (get (end chrom) pos)]
        (oc/load orid)))))


(defnp register-fragment
  "creates a fragment, registers and links it within the chromosome database
  performs a check before creation to avoid duplicates
  [data] hash-map of fragment properties"
  [{:keys [prev next chr p3 p5 dir]}]    ; dir in {:plus, :minus}, hash-maps {prev, next}
;  (println "register -> chr:" chr "p5:" p5 "\tp3:" p3 "\tnext:" next "  prev:" prev)
  (let [frag (or (select-fragment chr :p5 p5 dir)
                 (p :create-fragment (oc/with-tx
                   (oc/save! (og/vertex :fragment {:p3 p3, :p5 p5, :chr chr :dir dir})))))
        chr (select-chromosome chr)
        [l5 l3] (if (= dir :plus) [:p5frags :p3frags] [:mp5frags :mp3frags])]
    (if next
      (if-let [nfrag (select-fragment (:chr next) :p5 (:p5 next) (:dir next))]
        (link-frags frag nfrag)))
    (if prev
      (if-let [pfrag (select-fragment (:chr prev) :p3 (:p3 prev) (:dir prev))]
        (link-frags pfrag frag :count)))
    (link-with-chr frag l5 p5 chr)
    (link-with-chr frag l3 p3 chr)
     (assoc frag :next next :prev prev)))
; clj-orient bug: can't create line-break after call to :#rid


(defn remember-multistrand
  [frag]
  (if-not (= (:chr frag) (:chr (:next frag))) ; multistrand if elements on different strands
    (dosync
     (alter multis conj (:#rid frag)))) frag)


(defn remember-circular
  [frag]
  (if-let [next-p5 (:p5 (:next frag))]
    (let [is-upstream (if (= :plus (:dir frag)) < >)]
      (if ((every-pred true?)
               (= (:chr frag) (:chr (:next frag)))    ; circular if on same strand and
               (is-upstream next-p5 (:p5 frag)))      ; next frag starts upstream on same strand
        (dosync
         (alter circulars conj (:#rid frag)))))) frag)


;; Traversing the database requieres casting to and from record-ids
;; as different fetches of the same element may or may not be unique

(defn get-neighbours-
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
          (if-let [next-el (first (clojure.set/difference (get-neighbours- this) known))]
            (lazy-seq
             (cons this
                   (walk (conj known this) next-el)))))] (walk #{} start)))


(deftest structuretest
  (if (oc/db-exists? "memory:test")
    (do (println "\nDeleting database...")
        (oc/delete-db! (og/open-graph-db! "memory:test" "admin" "admin"))))
  (create-genome "memory:test" "/home/hermann/5_out.sam" {:head? #(complement (nil? %))
                                                          :header-parser #(do (identity %) {:id "2" :len 5000})
                                                          :data-parser identity
                                                          :add-all #(identity %)})
  (oc/with-db (og/open-graph-db! "memory:test" "admin" "admin")
    (is (nil? (select-chromosome "1")))
    (is (map? (add-chromosome "1" 4324)))
    (is (map? (select-chromosome "1")))
    (is (nil? (select-chromosome "error")))
    (let [f1 (register-fragment {:chr "1" :p5 1 :p3 11 :dir :plus :next {:chr "1" :p5 17 :dir :plus}})
          f2 (register-fragment {:chr "1" :p5 17 :p3 25 :dir :plus :prev {:chr "1" :p3 11 :dir :plus}})]
      (is (= 1 (count (og/get-links f1 f2)))))
    (is (nil? (select-fragment "1" :p5 17 :minus)))
    (is (map? (select-fragment "1" :p5 17 :plus)))))


(deftest sam-test
  (if (oc/db-exists? "memory:test")
    (do (println "\nDeleting database...")
        (oc/delete-db! (og/open-graph-db! "memory:test" "admin" "admin"))))
  (create-genome "memory:test" "test/data/5_out.sam"
                 {:head?         pressspan.saminput/header-line?
                  :header-parser pressspan.saminput/parse-header-line
                  :data-parser   pressspan.saminput/make-frag
                  :add-all (comp remember-circular remember-multistrand register-fragment)})
  (oc/with-db (og/open-graph-db! "memory:test" "admin" "admin")
    (is (= 5  (count (oq/native-query :chromosome {}))))
    (is (= 14 (count (oq/native-query :fragment {}))))
    (is (< 0 (og/count-edges)))
    (is (= 3 (count (get-subgraph (first @multis))))))
  (is (= 12  (count @multis))))

; (run-tests)

(if false
  (oc/with-db (og/open-graph-db! "memory:test" "admin" "admin")
    (for [frag (sort #(compare (:chr %1) (:chr %2)) (oq/native-query :fragment {}))]
      (println (format "%s:%d-%d" (:chr frag) (:p5 frag) (:p3 frag))))))
