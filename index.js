$.ajax({
    url: "https://quiz-app-367510.uc.r.appspot.com/quiz",
    type: 'GET'
}).then(function(r){
    console.log(r);
    document.body.appendChild(r);
});
