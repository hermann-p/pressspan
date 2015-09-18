(ns pressspan.graph
  (:require [clojure.data.int-map :as im]))


(use '[clojure.test]
     '[taoensso.timbre.profiling])


(defnp id-valid? [root id] (and (number? id) (>= id 0) (<= id (:nf @root))))


(defn- get-frag-id
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


(defn- create-chr
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


(defnp create-el
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


;(defn- fraginfo
;  [f]
;  ((juxt :chr :p5 :p3 :dir) f))
;(defn- same-frag?
;  [A B]
;  (zero? (compare (fraginfo A) (fraginfo B))))

(defn- is-link? [root up-id dn-id link-id]
  (let [link (get-in @root [:links link-id])]
    (and (= (:down link) dn-id)
         (= (:up link) up-id))))
(defnp get-link-id [root up-id dn-id]
  (assert (id-valid? root up-id) "Invalid upstream-fragment-id")
  (assert (id-valid? root dn-id) "Invalid downstream-fragment-id")
  (first (filter 
           (partial is-link? root up-id dn-id)
           (get-in @root [:frags up-id :down]))))
(defn- get-link [root up-id dn-id]
  (get-in @root [:links (get-link-id root up-id dn-id)]))


(defnp link-frags
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


(defnp remember-multistrand
  [root frag]
  (if-not (= (:chr frag) (:chr (:next frag))) ; multistrand if elements on different strands
    (dosync
     (commute update-in root [:multis] conj (:id frag)))) 
  frag)


(defnp remember-circular
  [root frag]
  (if-let [next-p5 (:p5 (:next frag))]
    (let [is-upstream? (if (= :plus (:dir frag)) < >)] ; define upstream-test
      (if ((every-pred true?)
               (= (:chr frag) (:chr (:next frag)))     ; circular if on same strand and
               (is-upstream? next-p5 (:p5 frag)))      ; next frag starts upstream on same strand
        (dosync
          (commute update-in root [:circulars] conj (:id frag))))))
  frag)


(defnp register-fragment
  [root {:keys[chr p3 p5 dir next prev] :as frag-data}]
  (let [frag (or (get-el root frag-data)
                 (create-el root frag-data))
        frag-id (:id frag)]
    (if prev
      (if-let [prev-id (:id (get-el root prev))]
        (link-frags root prev-id frag-id :count)))
    (assoc frag-data :id (:id frag))))


(defn- parse-header
  [root lines funs]; [{:keys [head? header-parser]}]]
  (doseq [line lines :while ((:head? funs) line)]
    (if-let [chromosome ((:header-parser funs) line)]
      (do
        (println "Chromosome:" chromosome)
        (create-chr root (:id chromosome) (:len chromosome)))
      (println line)))
  root)


(defnp parse-data
  [root lines funs];[{:keys [head? data-parser add-all]}]]
  (doseq [line (drop-while (:head? funs) lines)]
    ((:add-all funs) root ((:data-parser funs) line))))

(defn create-genome
  [filename funs]
  (println "Reading from file:" filename \newline)
  (let [genome (ref {:frags (im/int-map), :nf 0
                     :links (im/int-map), :nl 0
                     :circulars (im/int-set), :multis (im/int-set)})
        file (pressspan.io/lazy-file filename)]
    (assert (seq? file))
    (println "Processing header" \newline)
    (time (parse-header genome file funs))
    (println "Processing data...")
    (time (parse-data genome file funs))
    (println "File read.")
    @genome))


(deftest structuretest
  (let [A {:chr "1" :p5 1 :p3 177 :dir :plus}
        B {:chr "2" :p5 199 :p3 521 :dir :plus}
        C {:chr "2" :p5 1 :p3 127 :dir :plus}
        genome (ref {:frags (im/int-map), :nf 0, :links (im/int-map), :nl 0})]
    (create-chr genome "1" 1000)
    (create-chr genome "2" 1000)
    (create-el genome {:chr "1" :p5 200 :p3 300 :dir :plus})
    (is (map? (create-el genome A)))
    (is (map? (create-el genome B)))
    (is (map? (register-fragment genome (assoc C :prev B))))
    (is (map? (get-el genome A)))
    (is (map? (get-el genome B)))
    (is (nil? (get-el genome (assoc B :dir :minus))))
    (is (nil? (get-el genome (assoc A :chr "negative-test"))))
    (link-frags genome (:id (get-el genome A)) (:id (get-el genome B)))
    (link-frags genome (:id (get-el genome A)) (:id (get-el genome C)))))
    ;(println (clojure.string/join \newline @genome))))

(deftest samtest
  (let [filename "test/data/large.sam"
        funs {:head? pressspan.saminput/header-line?
              :header-parser pressspan.saminput/parse-header-line
              :data-parser pressspan.saminput/make-frag
              :add-all pressspan.graph/register-fragment}
        genome (create-genome filename funs)]
    (is (map? genome))
    (is (= 14 (:nf genome)))
    (is (= 10 (:nl genome)))))
 ;;   (println (clojure.string/join \newline genome))))