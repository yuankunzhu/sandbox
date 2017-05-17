var token, uid
var projectID = {}
var fileID = {}

jQuery.fn.extend({
    disable: function(state) {
        return this.each(function() {
            this.disabled = state;
        });
    }
});
$(document).ready(function() {
    $('#inputToken').bind('input propertychange', function() {
        $.getJSON('/api/v0/users', 
        {token: $('input[name="token"]').val()}, 
        function(data) {
            if(data.codes==200){
                token = $('input[name="token"]').val()
                uid = data.response.username
                $("#tokenBtn").removeClass('btn-danger').addClass('btn-success').text('valid');
                var description = "Authenticated user <b>" + data.response.first_name + "</b> from <b>" + data.response.affiliation + "</b>";
                $("#tokenDescription").html(description);
                $("#projectBtn").disable(false)
            } else {
                $("#tokenBtn").removeClass('btn-success').addClass('btn-danger').text('invalid');
                $('#tokenDescription').html("Invalid CAVATICA was provided, please input an correct one to move forward. You can find your authentication token from <a href='https://cavatica.sbgenomics.com/account/#developer'>Developer Dashboard</a>. See <a href='http://docs.cavatica.org/docs/the-api#section-authentication'>more details</a>.")
                $("#projectBtn").disable(true)
            }
        });
        return false;
    });
    $('#projectBtn').click(function() {
        $.getJSON('/api/v0/users/sbg-project', 
        {token: token}, 
        function(data) {
            var dataset = [];
            $.each(data.response, function(index, item) {
                var each = [];
                // project.push(item.id);
                projectID[item.name]=item.id
                each.push(item.name);
                dataset.push(each);
            });
            var table = $('#projects').DataTable( {
                destroy: true,
                pageLength: 5,
                bLengthChange: false,
                data: dataset,
                columns: [
                    //{ title: "ID" },
                    { title: "NAME" },
                    { title: "FILES" }
                ],
                columnDefs: [ {
                    targets: -1,
                    data: null,
                    defaultContent: "<button>file list</button>"
                }],
            });
            $('#projects tbody').on( 'click', 'button', function () {
                var pid = table.row( $(this).parents('tr') ).data();
                $.getJSON('/api/v0/users/sbg-project/files', 
                    {token: token, user: uid, project: projectID[pid[0]]}, 
                    function(data) {
                        var dataset = [];
                        $.each(data.response, function(index, item) {
                        var each = [];
                        fileID[item.name]=item.id
                        each.push(item.name);
                        dataset.push(each);
                    });
                    var table = $('#files').DataTable( {
                        destroy: true,
                        scrollY: 500,
                        sleect: true,
                        paging: false,
                        buttons:[
                            {text:'Select all',action: function(){table.rows().select}},
                            {text:'Select none',action: function(){table.rows().deselect}}
                        ],
                        data: dataset,
                        columns: [
                           { title: "NAME" }
                        ],
                        columnDefs: [ {
                            orderable: false,
                            className: 'select-checkbox',
                            targets: 0
                        }],
                        select:{style:'os',selector: 'td:not(:last-child)'}
                    });
                });
                // alert( "project id is " + projectID[data[0]] );
            });
        });
    });
    
});