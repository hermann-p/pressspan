(ns pressspan.saminput
  (:require [pressspan.io])
  (:use [clojure.string :only [split join]]
        clojure.test
;        taoensso.timbre.profiling))
        ))

(defn parse-header-line
  [line]
  (let [tokens (split line #"\s+")]
    (if (= "@SQ" (first tokens))
      (let [subtok (for [tok tokens] (split tok #":"))
            ids (for [el subtok] (first el))
            vals (for [el subtok] (last el))
            data (zipmap (map keyword ids) vals)]
        {:id (:SN data), :len (Integer. (:LN data))}))))


(defn header-line?
  [line]
  (= \@ (first line)))


(defn- extract-data
  "maps interesting fields of a tokenized alignment line from a .sam file"
  [tokens]
  (let [[_ flag rname p5 _ cigar _ _ _ _ _] tokens          ; destructure tokens
        sm-keys (re-seq #"X[XYQPUSTCV]:\w:[0-9A-Za-z.]+" (join " " (drop 11 tokens)))]  ; extract segemehl-info by regex
    {:chr rname :p5 (Integer. p5) :flag (Integer. flag) :cigar cigar :segemehl sm-keys})) ; return mapped results


(defn- extract-split-info
  "extracts split information from segemehl custom tags"
  [tokens]
  (let [tok (for [token tokens] (split token #":"))
        ids (for [el tok] (first el))
        vals (for [el tok] (last el))]
    (zipmap (map keyword ids) vals)))

(defn- cigar-ops [ops freqs]
  (reduce
     (fn [m tuple]
       (update m (keyword (first tuple)) (fn [k x] (+ (or k 0) (Integer. x))) (second tuple)))
     {}
     (map #(list %1 %2) ops freqs)))

(defn- un-cigar
  "calculates the original read length of an alignment with INDELs"
  [cigar]
  (let [freqs (re-seq #"[0-9]+" cigar)
        ops   (re-seq #"[MIDNHPX]" cigar)
        opvals (cigar-ops ops freqs)]
    (+ (:M opvals)
       (or (:D opvals) 0)
       (- (or (:I opvals) 0)))))



(defn make-frag
  "creates a read fragment node from a string of segemehl-created .sam data"
  [data-string]
  (let [tokens (split data-string #"\t+")
        data (extract-data tokens)
;        len (un-cigar (:cigar data))
        sm (extract-split-info (:segemehl data)) ; segemehl custom split data
        forward (= 0 (bit-and 16 (:flag data)))]  ; 0x10 => reverse
    
    {:dir (if forward :plus :minus)
     :chr (:chr data)
;         :p5 (Integer. (:XX sm))   ; segemehl reports positions in query, not 
;         :p3 (Integer. (:XY sm))}  ; on chromosome!

     :p5 (:p5 data)
     :p3 (+ (:p5 data) (un-cigar (:cigar data)) (- 1))
     
    	:prev (if-let [pchr (:XP sm)]
           {:chr pchr
            :p3 (Integer. (:XU sm))
            :dir (if (= "0" (:XS sm))
                   :minus
                   :plus)})

    	:next (if-let [nchr (:XC sm)]
        	{:chr nchr
         		:p5 (Integer. (:XV sm))
         		:dir (if (= "0" (:XT sm))
                :minus
                :plus)})
     }))

(deftest sam-processing-test
  (let [hdr-line "@SQ\tSN:1\tLN:3000"
        dta-line "Read1\t0\t1\t722\t255\t40M\t*\t0\t0\tACATCGAATCCACTACCAGTGAAACTGTGCCACAGCCACA\t*\tNM:i:0\tMD:Z:40\tNH:i:1\tXI:i:0\tXL:i:4\tXA:Z:\tQXX:i:41\tXY:i:80\tXQ:i:1\tXP:Z:1\tXU:i:40\tXS:i:64\tXC:Z:2\tXV:i:17\tXT:i:32
"]
    (is (= 10 (un-cigar "10M") (un-cigar "10M1I1D") (un-cigar "10D5M5I")))
    (is 33 (un-cigar "10M4D15M3I7M"))
    (is (= true (header-line? hdr-line)))
    (is (= false (header-line? dta-line)))
    (is (= {:id "1" :len 3000} (parse-header-line hdr-line)))
    (is (= {:dir :plus, :chr "1", :p5 722, :p3 761, :prev {:chr "1", :p3 40, :dir :plus}, :next {:chr "2", :p5 17, :dir :plus}} ; from manually created testdata
           (make-frag dta-line)))))

; (run-tests)
