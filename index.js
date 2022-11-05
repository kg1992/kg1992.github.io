const remoteUrl = "https://quiz-app-367510.uc.r.appspot.com/quote";
const localUrl = "http://localhost:8080/quote";

console.log("send request to localhost nodejs server")
$.ajax({
    url: localUrl,
    type: 'GET'
}).then(function(r){
    console.log(r);
    document.body.appendChild(r);
});

console.log("send request to remote google cloud app engine hosted nodejs server")
$.ajax({
    url: remoteUrl,
    type: 'GET'
}).then(function(r){
    console.log(r);
    document.body.appendChild(r);
});
