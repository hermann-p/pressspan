(ns pressspan.types)

(def flags {:PROCESSED 1 :CIRCULAR 2 :MULTISTRAND 4 :REVERSE 8})

(declare new-frag
         new-genome)

(definterface IFrag
  (setFlag [mode])
  (clrFlag [mode])
  (isSet [mode])
  (getUpLinks [])
  (getDnLinks [])
  (addUpLink [link])
  (addDnLink [link])
  (getDnCounts []))


(deftype Fragment [chr p5 p3
                   ^:volatile-mutable flags
                   ^:volatile-mutable up-links
                   ^:volatile-mutable dn-links
                   ^:volatile-mutable dn-cnts]
  IFrag
  (setFlag [this mode] (set! flags (bit-or flags mode)))
  (clrFlag [this mode] (set! flags (bit-xor flags mode))) ; because true ^ true = false
  (isSet [this mode] (> (bit-and flags mode) 0))
  (getUpLinks [this] up-links)
  (getDnLinks [this] dn-links)
  (getDnCounts [this] dn-cnts)
  (addUpLink [this link]
    (let [lnk-key (.p5 link)]
      (if up-links ; up-links nonempty
        (if-not (get up-links lnk-key)
          (set! up-links
                (assoc up-links lnk-key link)))
        (set! up-links ; up-links empty
              {lnk-key link}))))
  (addDnLink [this link]
    (let [link-key (.p5 link)]
      (if dn-links
        (if (get dn-links link-key)
          (set! dn-cnts ; dn links nonempty and this link registered
                (assoc dn-cnts link-key (inc (get dn-cnts link-key))))
          (do ; dn links nonempty,  but this link not registered
            (set! dn-links
                  (assoc dn-links link-key link))
            (set! dn-cnts
                  (assoc dn-cnts link-key 1))))
        (do  ; dn links empty
          (set! dn-links {link-key link})
          (set! dn-cnts {link-key 1}))))))




(definterface IChromosome
  (addFrag [p5 p3 upstr dnstr])
  (getFrag [p5])
  (getFrags []))

(deftype Chromosome [id length ^:volatile-mutable frags]
  IChromosome
  (addFrag [this p5 p3 upstr dnstr]
    (if-not (.getFrag this p5)
      (set! frags (assoc frags p5 (new-frag p5 p3))))
    (let [frag (.getFrag this p5)]
      (if upstr (.addUpLink frag upstr))
      (if dnstr (.addDnLink frag dnstr))))
      
  (getFrag [this p5] (get frags p5))
  (getFrags [this] frags))




(definterface IGenome
  (fromFile [file-name])
  (getChr [id])
  (addChr [id length])
  (addChrParallel [_ id length])
  (showChr []))

(deftype Genome [^:volatile-mutable chrs] ; chromosomes
  IGenome
  (getChr [this id] (get chrs id))
  (addChr [this id length] (set! chrs (assoc chrs id (Chromosome. id length {}))) this)
  (addChrParallel [this _ id length] (.addChr this id length) this)
  (showChr [this] chrs))


;; Quick constructor functions

(defn new-frag
  "Creates a new fragment, assigns an unset flag and no 3' or 5' links"
  [chr p5 p3]
  (Fragment. chr p5 p3 0 nil nil nil))

(defn new-genome
  "Creates an empty genome"
  []
  (Genome. {}))
