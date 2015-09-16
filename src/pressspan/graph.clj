(ns pressspan.graph
  (:require [clojure.data.int-map :as im]))


(use '[clojure.test]
     '[taoensso.timbre.profiling])


(def multis (ref #{}))
(def circulars (ref #{}))


(defn- get-el
	[root {:keys [chr p3 p5 dir]}]
  (if-not (get-in @root [chr])
    (println "No such chromosome:" chr))
	(let [[lnk pos] (if p5
						  [(if (= dir :plus) :p5frags :mp5frags) p5]
						  [(if (= dir :plus) :p3frags :mp3frags) p3])]
   (get-in @root [chr lnk pos])))


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


(defn- create-el
	[root {:keys [chr p3 p5 dir]}]
	(let [new-el {:chr chr :p3 p3 :p5 p5 :dir dir}
	      [l3 l5] (if (= dir :plus) [:p3frags :p5frags] [:mp3frags :mp5frags])]
   (dosync (alter root
                  #(-> %
                     (update-in [chr l3] into [[p3 new-el]])
                     (update-in [chr l5] into [[p5 new-el]]))))
   new-el))


(defn- fraginfo
  [f]
  ((juxt :chr :p5 :p3 :dir) f))
(defn- same-frag?
  [A B]
  (zero? (compare (fraginfo A) (fraginfo B))))


(defn- link-pos
  [up-el dn-el]
  (println "Finding linkpos..." \newline up-el \newline dn-el)
  (let [sf? (partial same-frag? dn-el)
        downlinks (map :down (:down up-el))]
    (printf "# of downlinks:" (count downlinks))
    (.indexOf downlinks (first (filter sf? downlinks)))))


(defn- link-frags
  [root up-el dn-el]
  (let [[up-chr up-p3 up-p5 up-dir] ((juxt :chr :p3 :p5 :dir) up-el)
        [up-l3 up-l5] (if (= :plus up-dir) [:p3frags :p5frags] [:mp3frags :mp5frags])
        [dn-chr dn-p3 dn-p5 dn-dir] ((juxt :chr :p3 :p5 :dir) dn-el)
        [dn-l3 dn-l5] (if (= :plus dn-dir) [:p3frags :p5frags] [:mp3frags :mp5frags])]
    (dosync
 ;     (if-let [link (get-link root up-el dn-el)]
 ;       ()
        ;; else create new link
        (let [link {:up up-el :depth 1 :down dn-el}
              up-el (update-in up-el [:down] conj link)
              dn-el (update-in dn-el [:up] conj link)]
          (alter root
                 #(-> %
                    (assoc-in [up-chr up-l3 up-p3] up-el)
                    (assoc-in [up-chr up-l5 up-p5] up-el)
                    (assoc-in [dn-chr dn-l3 dn-p3] dn-el)
                    (assoc-in [dn-chr dn-l5 dn-p5] dn-el)))))))


(defn- register-fragment
	[root {:keys [chr p3 p5 dir next prev] :as frag}]
	(let [el (or (get-el root frag) (create-el root frag))]
   (-> root
     #(if-let [next-frag (get-el % next)]
        (link-frags % frag next-frag))
     #(if-let [prev-frag (get-el % prev)]
        (link-frags % prev-frag frag :count)))))


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


(deftest structuretest
  (def genome (ref {}))
  (let [A {:chr "1" :p5 1 :p3 177 :dir :plus}
        B {:chr "2" :p5 199 :p3 521 :dir :plus}
        C {:chr "2" :p5 1 :p3 127 :dir :plus}]
    (create-chr genome "1" 1000)
    (create-chr genome "2" 1000)
    (create-el genome C)
    (create-el genome {:chr "1" :p5 200 :p3 300 :dir :plus})
    (is (map? (create-el genome A)))
    (is (map? (create-el genome B)))
    (is (map? (get-el genome A)))
    (is (map? (get-el genome B)))
    (is (nil? (get-el genome (assoc B :dir :minus))))
    (is (nil? (get-el genome {:chr "wrong one"})))
    (is (same-frag? A A))
    (is (not (same-frag? A B)))
    (link-frags genome A B)
    (link-frags genome A C)
    (println (clojure.string/join \newline @genome))
    (println (link-pos (get-el genome A) (get-el genome B)))))