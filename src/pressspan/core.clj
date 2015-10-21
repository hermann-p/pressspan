(ns pressspan.core
  (:import [java.io BufferedReader])
  (:require pressspan.io
            pressspan.saminput
            pressspan.graph)
   (:gen-class))

(require '[clojure.tools.cli :refer [parse-opts]])
;(require '[taoensso.timbre :as timbre
;           :refer (log  trace  debug  info  warn  error  fatal  report
;                        logf tracef debugf infof warnf errorf fatalf reportf
;                        spy get-env log-env)]
;         '[taoensso.timbre.profiling :as profiling
;           :refer (pspy pspy* profile defnp p p*)])
(use '[taoensso.timbre.profiling])


(def cli-options
  [["-m" "--multistrand" "Search for multistrand splits"]
   ["-c" "--circular" "Search for circular transcripts"]
   ["-i" "--in PATH_TO_FILE" "Input file, required"
    :default-desc "REQUIRED"]
   ["-o" "--out PATH_TO_FILE" "Output file name prefix"
    :default "pressspan_analysis"
    :default-desc "pressspan_analysis"]
   ["-r" "--report PATH_TO_FILE" "Create csv report and write it to FILE"
    :default "report.csv"
    :default-desc "report.csv"]
   ["-p" "--position CHR:xx [CHR:xx...]" "Output reads with fragments matching the  positions"
;    :default "0:0"
    :default-desc "N >= 1"]
   ["-R" "--range CHR:start:stop" "Output reads with fragments matching the range. Can only use -p or -R"]
   ["-d" "--drop DEPTH" "Drop isoform trees containing links with read depth lower than DEPTH"
;    :default 1
    :parse-fn #(Integer/parseInt %)
    :validate [#(< 0 %) "Must be a positive number"]]
   ["-t" "--trunc DEPTH" "Truncate isoform tree branches linked with read depth lower than DEPTH"
    :default 1
    :parse-fn #(Integer/parseInt %)
    :validate [pos? "Must be a positive number"]]
 ;  ["-D" "--database TYPE:NAME" "OrientDB database to use. Type in [local plocal remote]. Type remote needs running OrientDB server."
 ;   :default "memory:data"
 ;   :default-desc "memory"]
 ;  ["-U" "--user USERNAME" "Username for remote database"
 ;   :default "admin"]
 ;  ["-P" "--pass PASSWORD" "Password for remote database"
 ;   :default "admin"]
   ["-h" "--help"]])

(def functions
  {"sam" {:head? pressspan.saminput/header-line?
          :header-parser pressspan.saminput/parse-header-line
          :data-parser pressspan.saminput/make-frag
          :add-all pressspan.graph/register-fragment}})


(defn usage [summary]
  (println (clojure.string/join \newline [
   " "
   "pressspan 0.1"
   " "
   summary
   " "
   ])))

(defn error-text [errors summary]
  (usage summary)
  (println (clojure.string/join \newline errors)))


(defn exit-after [fun & stat]
  (try (fun) (catch Exception e))
  (System/exit (or stat 1)))



  
(defn -main
  [& args]
  (if-not [seq args] 
    (do
      (println (clojure.string/join \newline cli-options))
      (System/exit 1)))
  (let [{:keys [options arguments summary errors]} (parse-opts args cli-options)]
    (cond         ; check argument validity
;     (not (seq args)) (exit-after ((usage summary) summary 0)) ; duplicate...
     (:help options) (exit-after ((usage summary) summary 0))

     (< 0 (count errors))
     (exit-after (error-text errors summary))

     (nil? (:in options))
     (exit-after (error-text ["No input file --in given"] summary))

;     (every-pred (complement nil?) (:position options) (:range options))
                                        ;     (exit-after (error-text ["Can only select by --range or --position, else I'm confused, sorry"] summary))

     )

    (if-let [funs (get functions (last (clojure.string/split (:in options) #"\.")))]
      (let [input-funs [(if (:multistrand options) pressspan.graph/remember-multistrand)
                        (if (:circular options) pressspan.graph/remember-circular)
                        (:add-all funs)]
            funs (assoc funs :add-all input-funs)]
        (->
          (time (pressspan.graph/create-genome (:in options) funs))
          (pressspan.visualise/write-files :multis (:out options) (:trunc options))
          (pressspan.visualise/write-files :circulars (:out options))
          (pressspan.visualise/write-files :custom (:out options))))
      (println "No functions known to treat" 
               (last (clojure.string/split (:in options) #"\.")) ; get file extension
               "format. Why don't you create them?"))
    ))


;; Profiling test code
(if false
  (profile :info
           :Arithmetic
           (dotimes [n 1]
             (p :pressspan (-main "-i" "test/data/large.sam" "-m" "-c" )))))
