(ns pressspan.core
  (:import [java.io BufferedReader])
  (:require pressspan.io
            pressspan.saminput
            pressspan.graph
            pressspan.visualise
            pressspan.statistics
            [clojure.tools.cli :refer [parse-opts]])
  (:gen-class))

;; Structure for command line options
(def cli-options
  [["-m" "--multistrand" "Search for multistrand splits"]
   ["-c" "--circular" "Search for circular transcripts"]
   ["-i" "--in PATH_TO_FILE" "Input file, required"
    :default-desc "REQUIRED"]
   ["-o" "--out PATH_TO_FILE" "Output file name prefix"
    :default "pressspan_analysis"
    :default-desc "pressspan_analysis"]
   ["-S" "--stats" "Write event statistics CSV files"]
   ["-d" "--drop DEPTH" "Drop isoform trees containing links with read depth lower than DEPTH"
    :default 1
    :parse-fn #(Integer/parseInt %)
    :validate [#(< 0 %) "Must be a positive number"]]
;   ["-S" "--stats" "Write frequency of transsplicing events to file."]
   ["-s" "--bins N" "Define number of bins for frequency mapping"
    :default 10000
    :parse-fn #(Integer/parseInt %)
    :validate [pos? "Must be a positive number"]]
   ["-t" "--trunc DEPTH" "Truncate isoform tree branches linked with read depth lower than DEPTH"
    :default 1
    :parse-fn #(Integer/parseInt %)
    :validate [pos? "Must be a positive number"]]
   ["-h" "--help"]])

;; Maps file-type specific functions
(def functions
  {"sam" {:head? pressspan.saminput/header-line?
          :header-parser pressspan.saminput/parse-header-line
          :data-parser pressspan.saminput/make-frag
          :add-all pressspan.graph/register-fragment}})


(defn usage
  "Print a formatted help message"
  [summary]
  (println (clojure.string/join \newline [
                                          " "
                                          "pressspan 0.1"
                                          " "
                                          summary
                                          " "
                                          ])))

(defn error-text
  "Print error message and usage information"
  [errors summary]
  (usage summary)
  (println (clojure.string/join \newline errors)))


(defn exit-after
  "Execute a function, then exit cleanly. Called when errors occured."
  [fun & stat]
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
      (:help options)
      (exit-after ((usage summary) summary 0))

      (< 0 (count errors))
      (exit-after (error-text errors summary))

      (nil? (:in options))
      (exit-after (error-text ["No input file --in given"] summary)))
    
    (if-let [funs (get functions (last (clojure.string/split (:in options) #"\.")))]
      (let [input-funs [(if (:multistrand options) pressspan.graph/remember-multistrand)
                        (if (:circular options) pressspan.graph/remember-circular)
                        (:add-all funs)]
            drop-graphs (pressspan.visualise/drop-filter (:drop options))
            funs (assoc funs :add-all input-funs)]
        (->
         (time (pressspan.graph/create-genome (:in options) funs (:bins options)))
         (pressspan.visualise/write-files :multis (:out options) (:trunc options) [drop-graphs])
         (pressspan.visualise/write-files :circulars (:out options) (:trunc options) [drop-graphs])
         (#(when (:stats options)
             (pressspan.statistics/write-stat-file % (str (:out options) "/multistrand.csv") :multis)
             (pressspan.statistics/write-stat-file % (str (:out options) "/circulars.csv") :circulars)))))
      (println "No functions known to treat" 
               (last (clojure.string/split (:in options) #"\.")) ; get file extension
               "format. Why don't you create them?"))
    (shutdown-agents)))
