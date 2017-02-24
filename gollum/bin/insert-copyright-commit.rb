require 'rubygems'
require 'gollum-lib'

wiki = Gollum::Wiki.new('pinboard-bookmarks-to-chrome.wiki')

def create_commit(page)
  { :message => "added copyright to #{page}",
    :name => 'Grant Winney',
    :email => 'user@email.com' }
end

copyright = '<p><em>Copyright 2017 - Grant Winney - <a href="https://opensource.org/licenses/MIT">MIT License</a></em></p>'

Dir.foreach('pinboard-bookmarks-to-chrome.wiki') do |item|
  ext = File.extname(item)
  next if ext != '.md'
  basename = File.basename(item, ext)
  page = wiki.page(basename)
  if (!page.raw_data.start_with?(copyright))
    wiki.update_page(page,
                     page.name,
                     page.format,
                     "#{copyright}\n\n#{page.raw_data}",
                     create_commit(page.name))
  end
end