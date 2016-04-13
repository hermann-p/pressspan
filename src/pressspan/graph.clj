(ns pressspan.graph
  (:require [clojure.data.int-map :as im]
            [pressspan.statistics :as stats]
            [clojure.core.async
             :as a
             :refer [>! <! >!! <!! go chan close! thread
                     alts! alts!! timeout go-loop put!]]
            pressspan.io)
  (:use [clojure.set :only [difference union]]
        clojure.test
        ))

(def frag-chan (chan 10))
(def mult-chan (chan 10))
(def circ-chan (chan 10))
(declare mult-freqs
         circ-freqs)

(def parsing-done (atom 0))

(defn id-valid? [root id] (and (number? id) (>= id 0) (<= id (:nf @root))))


(defn get-frag-id
	[root {:keys [chr p3 p5 dir] :as query}]
	(let [[lnk pos] (if p5
                   [(if (= dir :plus) :p5frags :mp5frags) p5]
                   [(if (= dir :plus) :p3frags :mp3frags) p3])]
   (get-in @root [chr lnk pos])))


(defmulti get-el (fn [_ el] (map? el)))

(defmethod get-el true [root {:keys [chr p3 p5 dir id] :as query}]
  (if-not (get-in @root [chr])
    (println "No such chromosome:" chr)
    (if-let [id (or id (get-frag-id root query))]
      (get-el root id))))

(defmethod get-el false [root id]
  (assert (id-valid? root id) "Invalid node id")
  (get-in @root [:frags id]))


(defn create-chr
  [root id len]
  (dosync 
    (alter root
           #(-> %
              (assoc-in [id] {:len len})
              (update-in [id] into
                          [[:p3frags (im/int-map)]
                           [:p5frags (im/int-map)]
                           [:mp3frags (im/int-map)]
                           [:mp5frags (im/int-map)]])))))


(defn create-el
	[root {:keys [chr p3 p5 dir]}]
    (dosync 	
      (let [[l3 l5] (if (= dir :plus) [:p3frags :p5frags] [:mp3frags :mp5frags])
           id (:nf @root)
           new-el {:chr chr :p3 p3 :p5 p5 :dir dir :id id}]
        (alter root
               #(-> %
                  (update-in [chr l3] into [[p3 id]])
                  (update-in [chr l5] into [[p5 id]])
                  (update :frags into [[id new-el]])
                  (update :nf inc)))
        new-el)))


(defn is-link? [root up-id dn-id link-id]
  (let [link (get-in @root [:links link-id])]
    (and (= (:down link) dn-id)
         (= (:up link) up-id))))
(defn get-link-id [root up-id dn-id]
  (assert (id-valid? root up-id) "Invalid upstream-fragment-id")
  (assert (id-valid? root dn-id) "Invalid downstream-fragment-id")
  (first (filter 
           (partial is-link? root up-id dn-id)
           (get-in @root [:frags up-id :down]))))
(defn get-link [root up-id dn-id]
  (get-in @root [:links (get-link-id root up-id dn-id)]))


(defn link-frags
  [root up-id dn-id & count]
  (assert (id-valid? root up-id) "Invalid upstream-fragment-id")
  (assert (id-valid? root dn-id) "Invalid downstream-fragment-id")
  (if-let [link-id (get-link-id root up-id dn-id)]
      (if count (dosync (commute root update-in [:links link-id :depth] inc)))
    ;; else create new link
    (dosync
      (let [link-id (get-in @root [:nl])
            link {:down dn-id :depth (if count 1 0) :id link-id :up up-id}]
        (if-not (get-in @root [:frags up-id :down])
          (alter root assoc-in [:frags up-id :down] (im/int-set)))
        (if-not (get-in @root [:frags dn-id :up])
          (alter root assoc-in [:frags dn-id :up] (im/int-set)))
        (alter root
               #(-> %
                  (update-in [:nl] inc)
                  (update-in [:links] into [[link-id link]])
                  (update-in [:frags up-id :down] into [link-id])
                  (update-in [:frags dn-id :up] into [link-id])))))))


(defn remember-multistrand
  [genome frag]
  (when-let [next (:next frag)]
    (when-not (and
             (= (:chr frag) (:chr next))
             (= (:dir frag) (:dir next))); multistrand if elements on different strands
      (dosync
       (commute genome update-in [:multis] conj (:id frag)))
      (go (>! mult-chan [(:chr frag) (:p3 frag) (:dir frag)])
          (>! mult-chan [(:chr next) (:p5 next) (:dir next)]))))
  frag)


(declare circ-freqs)
(defn remember-circular
  [genome frag]
  (if-let [next-p5 (:p5 (:next frag))]
    (let [is-upstream? (if (= :plus (:dir frag)) < >)] ; define upstream-test
      (when ((every-pred true?)
      				 (= (:dir frag) (:dir (:next frag)))
               (= (:chr frag) (:chr (:next frag)))     ; circular if on same strand andh
               (is-upstream? next-p5 (:p5 frag)))      ; next frag starts upstream on same strand
        (dosync
         (commute genome update-in [:circulars] conj (:id frag)))
        (go (>! circ-chan [(:chr frag) (:p3 frag) (:dir frag)])
            (>! circ-chan [(:chr next) (:p5 next) (:dir next)])))))
  frag)


(defn register-fragment
  [root {:keys[chr p3 p5 dir next prev] :as frag-data}]
  (let [frag (or (get-el root frag-data)
                 (create-el root frag-data))
        frag-id (:id frag)]
    (if prev
      (if-let [prev-id (:id (get-el root prev))]
        (link-frags root prev-id frag-id :count)))
    (assoc frag-data :id (:id frag))))


(defn parse-header
  [root lines funs]; [{:keys [head? header-parser]}]]
  (doseq [line lines :while ((:head? funs) line)]
    (if-let [chromosome ((:header-parser funs) line)]
      (do
        (println "Chromosome:" chromosome)
        (create-chr root (:id chromosome) (:len chromosome)))
      (println line)))
  root)


(defn parse-data
  [root lines funs];[{:keys [head? data-parser add-all]}]]
  (let [add-frag (filter identity (:add-all funs))
        add-frag (mapv #(partial % root) add-frag)
        add-frag (reduce comp add-frag)]
    ;; explicitly specify tree-building thread
    (thread ;; extra thread because called for every frag
      (loop []
        (let [[frag _] (alts!! [frag-chan])]
          (if-not (= frag :done)
            (do (add-frag frag)
                (recur))
            (swap! parsing-done inc)))))

    ;; explicitly specify event-counting threads
    (go-loop []
      (let [[data _] (alts! [mult-chan])]
        (if-not (= data :done)
          (do (let [[chr pos dir] data]
                (stats/increase root mult-freqs chr pos dir))
              (recur))
          (swap! parsing-done inc))))
    (go-loop []
      (let [[data _] (alts! [circ-chan])]
        (if-not (= data :done)
          (do (let [[chr pos dir] data]
                (stats/increase root circ-freqs chr pos dir))
              (recur))
          (swap! parsing-done inc))))

    ;; feed parsed fragments to fragment channel
    (doseq [line (drop-while (:head? funs) lines)]
      (>!! frag-chan ((:data-parser funs) line)))
    ;; send threads ending signals
    (>!! frag-chan :done)
    (>!! mult-chan :done)
    (>!! circ-chan :done)
    (while (< @parsing-done 3) (Thread/sleep 10))))

(defn create-genome
  [filename funs & bins]
  (println "Reading from file:" filename \newline)
  (let [genome (ref {:frags (im/int-map), :nf 0
                     :links (im/int-map), :nl 0
                     :circulars (im/int-set), :multis (im/int-set)})
        file (pressspan.io/lazy-file filename)]
    (assert (seq? file))
    (println "Processing header" \newline)
    (parse-header genome file funs) ;; add header information to genome
    ;; Create new parallelly accessible data...
    (def circ-freqs (ref (stats/empty-stats @genome (or bins 10000))))
    (def mult-freqs (ref (stats/empty-stats @genome (or bins 10000))))
    (println "Processing data...")
    (time (parse-data genome file funs))
    (println "File read.")
    (-> @genome
        (assoc :mstat @mult-freqs
               :cstat @circ-freqs))))


;;;#############################################################################
;;; Postprocessing functions work on maps, NOT ref's
;;;#############################################################################

(defn get-both [graph trunc lnk-id]
  (let [el (get-in graph [:links lnk-id])]
 ;   (println "el:" el "minimum link depth:" trunc "pass:" (<= trunc (:depth el)))
    (if (<= trunc (:depth el))
      ((juxt :up :down) el))))

(defn get-adjacent 
  ([graph el-id] (get-adjacent graph el-id 1))
  ([graph el-id trunc]
;    (println "get-adjacent graph" el-id "trunc:" trunc)
    (let [el (get-in graph [:frags el-id])
        link-ids (apply concat [(:up el) (:down el)])
        links (map (partial get-both graph trunc) link-ids)
        links (set (apply concat links))]
      (difference links #{el-id}))))


(defn prune [root node-id min-depth]
  (let [node (get-in root [:frags node-id])
        up-links (map #(get-in root [:links %]) (:up node))
        dn-links (map #(get-in root [:links %]) (:down node))
        truncate (fn [links] (set (filter #(<= min-depth (:depth %)) links)))]
 ;   (println "Node:" node "===> "
    (->
      node
      (assoc :up (truncate up-links))
      (assoc :down (truncate dn-links)))
    (->
      node
      (assoc :up (truncate up-links))
      (assoc :down (truncate dn-links)))))


(defn get-subgraph 
  ([graph start] (get-subgraph graph start 1))
  ([graph start min-depth]
    (let [walk
        (fn walk [known unvisited this]
;          (println "walk\nknown:" known "unvisited:" unvisited "this:" this)
          (let [unvisited (union (get-adjacent graph this min-depth) unvisited)
                unvisited (difference unvisited known)
                unvisited (difference unvisited #{this})]
            (if-let [next-el (first (difference unvisited known))]
              (lazy-seq
                (cons (prune graph this min-depth)
                      (walk (conj known this) (difference unvisited #{this}) next-el)))
              (lazy-seq [(prune graph this min-depth)]))))]
      (walk #{} #{} start))))


(defn all-subgraphs 
  ([graph indices] (all-subgraphs graph indices 1))
  ([graph indices trunc]
    (let [walk
        (fn walk [unvisited]
          (if-let [start (first unvisited)]
            (let [sg (get-subgraph graph start trunc)]
              (if (seq sg)
                (lazy-seq
                  (cons sg
                        (walk (difference unvisited (set (map :id sg))))))))))]
      (walk indices))))
            
