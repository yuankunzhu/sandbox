require 'rubygems'
require 'gollum-lib'

wiki = Gollum::Wiki.new('pinboard-bookmarks-to-chrome.wiki')
html = '<h1>Wiki Contents</h1>'
html << '<style>p { font-size:larger; font-weight:bold }</style>'

Dir.foreach('pinboard-bookmarks-to-chrome.wiki') do |item|
  ext = File.extname(item)
  next if ext != '.md'
  basename = File.basename(item, ext)
  page = wiki.page(basename)
  html << "<hr><p><a href=\"https://github.com/grantwinney/pinboard-bookmarks-to-chrome/wiki/#{basename}\">#{page.name}</a></p>"
  toc = page.toc_data
  next if toc.nil?
  html << "#{toc}"
end

File.write('outline.html', html)