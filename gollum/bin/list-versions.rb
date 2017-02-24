require 'rubygems'
require 'gollum-lib'

wiki = Gollum::Wiki.new('pinboard-bookmarks-to-chrome.wiki')
Dir.foreach('pinboard-bookmarks-to-chrome.wiki') do |item|
  ext = File.extname(item)
  next if ext != '.md'
  page = wiki.page(File.basename(item, ext))
  puts "#{item}: #{page.version.id}\n"
end