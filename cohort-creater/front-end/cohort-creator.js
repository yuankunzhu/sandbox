var projecid = {};
var sampleid = {};
var filelist = [];
var fulllist = [];

jQuery.fn.extend({
    disable: function(state) {
        return this.each(function() {
            this.disabled = state;
        });
    }
});
function getUserInfo(){
    $.getJSON('./data/user.json', function(data){
        document.getElementById("user_name").outerHTML=data.name;
        document.getElementById("user_email").outerHTML=data.email;
        document.getElementById("user_affiliation").outerHTML=data.affiliation;
    }); 
};
function validateForm(str) {
    var allLetters = /^[a-zA-Z0-9\-\_]+$/;
    if (allLetters.test(str)) {
        $('#runBtn').text("Run Cohort!")
        $('#runBtn').disable(false)
    }
    else {
        $('#runBtn').text("Invalid Naming Letter")
        $('#runBtn').disable(true)
    }
};
$(document).ready(function() {
    $('#projectBtn').click(function() {
        $("#debug").text("[READ]: data/project.json");
        filelist = [];
        $.getJSON('./data/project.json', function(data){
            var dataset = [];
            $.each(data, function(index, item) {
                var each = [];
                projecid[item.name]=item.id
                var dis = 'disabled>No Access';
                if(item.exe){
                    dis = '>File List';
                }
                var btn = '<button class="btn btn-xs" '+dis+'</button>';
                each.push(item.name);
                each.push(btn);
                dataset.push(each);
            });
            var table = $('#projects').DataTable( {
                bFilter: false, bInfo: false,
                destroy: true,
                paging: false,
                bLengthChange: false,
                data: dataset,
                columns: [
                    { title: "Name" },
                    { title: "Files" }
                ]
            });
            $('#projects tbody').on( 'click', 'button', function () {
                var tablerow = table.row( $(this).parents('tr') ).data();
                filelist = [];
                $("#debug").html("[RETURN project-id]: "+projecid[tablerow[0]]+"<br>[READ] data/files.json")
                $("#dropBtn").text("Pipeline");
                $("#cohortName").disable(true);
                $("#dropBtn").disable(true);
                $("#runBtn").disable(true);
                
                $("#projectBtn").text(tablerow[0])
                $.getJSON('./data/file.json', function(data){
                    var dataset = [];
                    $.each(data, function(index, item) {
                        var each = [];
                        var dis = '<i class="fa fa-ban" aria-hidden="true"></i>'
                        if(item.metadata){
                            dis = '<input type="checkbox">';
                        }
                        fulllist.push(item.id);
                        sampleid[item.name]=item.id;
                        each.push(item.name);
                        each.push(item.type);
                        each.push(item.metadata);
                        each.push(dis);
                        dataset.push(each);
                    });
                    var table = $('#files').DataTable( {
                        bInfo: false,
                        destroy: true,
                        paging: false,
                        bLengthChange: false,
                        data: dataset,
                        columnDefs: [{
                            className: "dt-center", 
                            targets: 3,
                            orderable: false
                        }],
                        columns: [
                            { title: "Name" },
                            { title: "Data Type" },
                            { title: "Metadata" },
                            { title: '<input type="checkbox" id="check-all">' }
                        ]
                    });
                    $('#check-all').on('click', function(){
                        var rows = table.rows({ 'search': 'applied' }).nodes();
                        $('input[type="checkbox"]', rows).prop('checked', this.checked);
                        var checkCount = 0;
                        table.$('input[type="checkbox"]').each(function(){
                            filelist = fulllist;
                            if(this.checked){
                                checkCount+=1;
                            }
                        });
                        if (checkCount==0){
                            $("#cohortName").disable(true); 
                            $("#dropBtn").disable(true);
                            $("#runBtn").disable(true);
                        }else{
                            $("#cohortName").disable(false);
                            $("#dropBtn").disable(false);
                        }
                    });
                    $('#files tbody').on('change', 'input[type="checkbox"]', function(){
                        var checkCount = 0;
                        var tablerow = table.row( $(this).parents('tr') ).data();
                        filelist = [];
                        table.$('input[type="checkbox"]').each(function(){
                            if(this.checked){
                                checkCount+=1;
                                filelist.push(sampleid[tablerow[0]])
                            }
                        });
                        if (checkCount==0){
                            $("#cohortName").disable(true);
                            $("#dropBtn").disable(true);
                            $("#runBtn").disable(true);
                        }else{
                            $("#cohortName").disable(false);
                            $("#dropBtn").disable(false);
                        }
                        if(!this.checked){
                            // $("#dropBtn").disable(true);
                            var el = $('#check-all').get(0);
                            // If "Select all" control is checked and has 'indeterminate' property
                            if(el && el.checked && ('indeterminate' in el)){
                                // Set visual state of "Select all" control 
                                // as 'indeterminate'
                                el.indeterminate = true;
                            }
                        }
                        else{
                            $("#dropBtn").disable(false);
                        }
                    });
                });
            });
        });
    });
    $('#cohortName').bind('input propertychange', function() {
        validateForm($('input[name="cName"]').val())
    });
    $('#wgs').click(function() {
        $("#dropBtn").text("WGS")
    });
    $('#wes').click(function() {
        $("#dropBtn").text("WES")
    });
    $('#rna').click(function() {
        $("#dropBtn").text("RNA-Seq")
    });
    $('#runBtn').click(function() {
        $("#debug").html("[RETUN]: -Pipeline-: "+$("#dropBtn").text()+"<br>-cohort name-: "+$("#cohortName").val()+"<br>-file list-: "+filelist);
    });
});
