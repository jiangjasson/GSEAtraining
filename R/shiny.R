

#' Interactive Table for AnnotationHub
#' 
#' @import DBI
#' @import RSQLite
#' @export
ah_shiny = function() {

    check_pkg("AnnotationHub", bioc = TRUE)
    check_pkg("shiny")
    check_pkg("DT")
    check_pkg("htmltools")

    ah = AnnotationHub::AnnotationHub()
    sqlite = dbDriver("SQLite")
    db = dbConnect(sqlite, ah@.db_path)

    tb = dbGetQuery(db, "SELECT r.ah_id, r.title, r.species, r.taxonomyid, r.genome, r.description, p.rdataclass, r.rdatadateadded from resources r inner join rdatapaths p on r.id = p.id")

    dbDisconnect(db)

    ui = shiny::fluidPage(
        title = 'AnnotationHub',
        htmltools::h2("Interactive Table for AnnotationHub"),
        htmltools::hr(),
        DT::dataTableOutput('tbl')
    )

    server = function(input, output, session) {
        output$tbl = DT::renderDataTable(tb, server = TRUE, filter = "top")
    }

    print(shiny::shinyApp(ui, server))
}

