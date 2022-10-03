(defproject bioit-exam "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "EPL-2.0 OR GPL-2.0-or-later WITH Classpath-exception-2.0"
            :url "https://www.eclipse.org/legal/epl-2.0/"}
  :dependencies [[org.clojure/algo.generic "0.1.3"]
                 [org.clojure/clojure "1.11.1"]
                 [org.clojure/core.memoize "1.0.257"]]

  :profiles {:dev {:dependencies [[org.clojure/test.check "0.10.0"]]}}

  :repl-options {:init-ns bioit-exam.core})
