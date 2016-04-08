(defproject pressspan "0.1.2"
  :description "A segemehl reassembly and visualisation tool"
  :url "https://bitbucket.org/waechtertroll/pressspan"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.7.0"]
                 [org.clojure/data.int-map "0.2.0"]
                 [org.clojure/tools.cli "0.3.2"]]
  :main ^:skip-aot pressspan.core
  :target-path "target/%s"
  :profiles {:uberjar {:aot :all}})
