class HtmlTabWriter:
    
    def __init__(self, outDir):
        self.outDir = outDir
        self.html_file = f'{self.outDir}/index.html'
        self.fo_html = open(self.html_file, 'w')
        self.abundacy_tab = f'{self.outDir}/clust_abund.tab'
        self.annotation_dir = f'{self.outDir}/diamond_annotation'
        self.main()
    
    def _get_head(self):
        return """<html>
<head>
<style>
.container {
  padding: 2rem 0rem;
}

h4 {
  margin: 2rem 0rem 1rem;
}

.table-image {
  td, th {
    vertical-align: middle;
  }
}
.is-breakable {
  word-break: break-word;
}
table {
    table-layout: fixed;
    width:100%
}
td {
    word-wrap:break-word;
}

/*ZOOM MAPY*/
.zoom {

    transition: transform .08s;
  }
  
  .zoom:hover {
    -ms-transform: scale(1.5); /* IE 9 */
    -webkit-transform: scale(1.5); /* Safari 3-8 */
    transform: scale(5.5); 
  }
</style>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@4.3.1/dist/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
<link rel="stylesheet" href="//netdna.bootstrapcdn.com/bootstrap/3.0.0/css/bootstrap-glyphicons.css">
<script src="https://kit.fontawesome.com/1a74312711.js" crossorigin="anonymous"></script>


<script>
function SelfCopy(copyText)
  {
      navigator.clipboard.writeText(copyText);
      alert("Sequence copied!");
  }
</script>

<title>nanoTRF report</title>
</head>
<body>
<div class="container">
  <div class="row">
    <div class="col-12">
		<table class="table table-image">
		  <thead>
		    <tr>
		      <th scope="col">Number</th>
		      <th scope="col">Cluster Name</th>
              <th scope="col"> CLuster </th>
              <th scope="col">Genome abundancy, %</th>
              <th scope="col">Genome abundancy corrected, %</th>
              <th scope="col"> Per Read Coverage By Cluster, % </th>
              <th scope="col"> Per Read Coverage By Cluster Pie, % </th>
		      <th scope="col">Min monomer size assembled by Cap3</th>
		      <th scope="col">Max monomer size assembled by Cap3</th>
              <th scope="col">Subrepeats length</th>
              <th scope="col">Annotation</th>
              <th scope="col">Per read annotation</th>
		    </tr>
		  </thead>
		  <tbody>

"""
    
    def _get_footer(self):
        return 		 """ </tbody>
		</table>   
    </div>
  </div>
</div>
  
</body>
</html>
"""
        
    def main(self):
        
        self.fo_html.write(self._get_head())
        
        with open(self.abundacy_tab) as inFile:
            for i, lines in enumerate(inFile):
                sp = lines.rstrip().split('\t')
                if i!=0 and 'artef' not in sp[0]:
                    cluster,minContig,maxContig,abundancy,abundancy_corr, Contig1_seq = lines.rstrip().split('\t')[:6]
                    cluster_png = f'clusters/{cluster}.png'
                    cluster_png_hist = f'clusters/{cluster}.hist.png'
                    cluster_png_pie = f'clusters/{cluster}.pie.png'
                    cluster_png_annotation = f'diamond_annotation/{cluster}.dom.png'
                    if sp[-3] != 'nf':
                        trf_seq, trf_len, annotation = lines.rstrip().split('\t')[6:]
                        
                        ## take one hit only for html
                        if ";" in annotation:
                            annotation = annotation.split(';')[0]
                        min_len = int(trf_len.split(' - ')[0])
                        for seq in trf_seq.split(';'):
                            if len(seq.split(':')[1]) == min_len:
                                Contig1_seq = seq
                        Contig1_seq = f"{min_len}bp:{Contig1_seq}"

                    tr = f"""<tr>
                                  <th scope="row">{i}</th>
                                  <td>{cluster}</td>
                                  <td class="w-25"">
                                      <div class='zoom'>
                                          <img src="{cluster_png}" class="img-fluid img-thumbnail" alt="{cluster}">
                                      </div>
                                  </td>
                                  <td> 
                                          {abundancy}
                                  </td>
                                  <td> 
                                          {abundancy_corr}
                                  </td>
                                  <td class="w-25"">
                                      <div class='zoom'>
                                          <img src="{cluster_png_hist}" class="img-fluid img-thumbnail" alt="{cluster}">
                                      </div>
                                  </td>
                                  <td class="w-25"">
                                   <div class='zoom'>
                                      <img src="{cluster_png_pie}" class="img-fluid img-thumbnail" alt="{cluster}">
                                  </div>
                                  </td>
                                  <td >{minContig}</td>
                                  <td>{maxContig}</td>
                                  <td> 
                                          {trf_len}
                                  </td>
                                    <td>{annotation}</td>
                                  <td class="w-25"">
                                   <div class='zoom'>
                                      <img src="{cluster_png_annotation}" class="img-fluid img-thumbnail" alt="{cluster}">
                                  </div>
                                  </td>


                </tr>"""
                    print(cluster_png)
                    self.fo_html.write(tr)
        
        self.fo_html.write(self._get_footer())
        
        self.fo_html.close()
    