
//LLIBRERIA
const KB = [
  {
    id: 1,
    keys: ["oïdi","erysiphe","necator","cendrosa","fong","malaltia","patogen","mildiu","espora","fungi","cep","vinya"],
    title: "Erysiphe necator (Oïdi / Cendrosa)",
    body: "L'<b>Erysiphe necator</b>, popularment conegut com a <em class='term'>cendrosa</em> o oïdi, és el patogen fúngic més persistent de la viticultura mediterrània. A diferència del míldiu, no depèn de la pluja: necessita únicament humitat ambiental i temperatures moderades (15–28°C) per propagar-se. L'oïdi pot destruir el 100% d'una collita en atacs severs i contamina el raïm amb compostos volàtils que inhabiliten la producció de vi de qualitat.",
    highlight: "Esporula sense pluja, forma colònies blanques polvoroses sobre els teixits verds, i sobreviu l'hivern en forma de cleistotècies als sarments.",
    suggestions: ["compostos volàtils fong","resistència fungicides sofre","canvi climàtic temperatura"],
    confidence: 0.97
  },
  {
    id: 2,
    keys: ["biosensor","sensor","detecció","detectar","dispositiu","nas","electrònic","com","funciona","funcionament","mecanisme","principi","explica","sistema","tecnologia","opera","treballa"],
    title: "El Biosensor: Com Funciona",
    body: "El nostre biosensor actua com un <b>nas electrònic biològic</b>. Integra proteïnes d'unió a odors (<em class='term'>OBPs</em>) en un transductor electroquímic. Quan l'oïdi comença a infectar el cep, allibera compostos volàtils específics que s'uneixen a les OBPs. Aquesta unió provoca un canvi conformacional que el sensor tradueix en un <b>senyal elèctric mesurable</b>, enviat en temps real al sistema de gestió del viticultor.",
    highlight: "El sensor detecta el fong ABANS que el dany sigui visible a ull nu — és detecció precoç, no diagnòstic tardà.",
    suggestions: ["proteïnes OBP receptor","compostos volàtils fong","senyal elèctric transducció"],
    confidence: 0.99
  },
  {
    id: 3,
    keys: ["obp","proteïna","proteines","olfactory","binding","protein","insecte","antena","receptor","unió","conformació","lligand","recombinant","disseny"],
    title: "Proteïnes OBP (Olfactory Binding Proteins)",
    body: "Les <b>Olfactory Binding Proteins (OBPs)</b> són proteïnes presents en les antenes d'insectes que capturen i transporten molècules odoroses fins als receptors olfactoris. El nostre equip ha dissenyat OBPs <em class='term'>recombinants</em> amb especificitat cap als compostos volàtils característics d'<em>Erysiphe necator</em>. En unir-se al compost, l'OBP canvia de conformació (plegament 3D), generant un senyal detectable.",
    highlight: "Clau del projecte: mimetitzem la natura. Els insectes porten milions d'anys detectant olors amb altíssima sensibilitat i selectivitat.",
    suggestions: ["compostos volàtils fong","com funciona el biosensor","enginyeria dbtl"],
    confidence: 0.98
  },
  {
    id: 4,
    keys: ["volàtil","volatils","compostos","aroma","olor","molècula","cov","hexanol","hexanal","terpè","cariofile","ppb","química","empremta","quimica"],
    title: "Compostos Volàtils d'Erysiphe necator",
    body: "Quan <em>Erysiphe necator</em> colonitza el cep, emet un perfil específic de <b>compostos orgànics volàtils (COVs)</b>: principalment <em class='term'>1-hexanol</em>, <em class='term'>hexanal</em>, i certs terpens fúngics com el β-cariofilè. Aquests compostos diferen dels emesos per la planta sana, cosa que permet discriminar la infecció activa de l'absència de patogen.",
    highlight: "La 'empremta química' del fong és única i detectable a concentracions de parts per bilió (ppb), molt abans que apareguin símptomes visuals.",
    suggestions: ["proteïnes OBP receptor","com funciona el biosensor","oïdi erysiphe cendrosa"],
    confidence: 0.96
  },
  {
    id: 5,
    keys: ["resistència","fungicida","sofre","tractament","control","quimic","pesticides","ibe","sdhi","fitotoxicitat","ecològic","ecologia"],
    title: "Resistència als Fungicides i Limitacions del Control Actual",
    body: "El control convencional de l'oïdi enfronta dos problemes majors: (1) <b>resistències genètiques</b> del fong a famílies d'IBE (inhibidors de biosíntesi d'ergosterol) i d'SDHI, reduint l'eficàcia dels fungicides sistèmics. (2) L'<b>agricultura ecològica catalana</b> depèn del sofre, però les onades de calor (&gt;35°C) provoquen fitotoxicitat i danya el cep, creant finestres d'inactivitat forçada precisament en moments de risc alt.",
    highlight: "Resultat: els viticultors han d'aplicar tractaments preventius sense saber si el fong és present, disparant els costos de producció.",
    suggestions: ["impacte econòmic collita vi","com funciona el biosensor","canvi climàtic temperatura"],
    confidence: 0.95
  },
  {
    id: 6,
    keys: ["econòmic","economia","cost","pèrdua","collita","producció","benefici","marge","vi","cava","viticultor","rendibilitat","diners","preu"],
    title: "Impacte Econòmic a la Viticultura Catalana",
    body: "L'oïdi és la malaltia que genera <b>majors costos i incertesa</b> al sector vitivinícola català. En atacs severs, la pèrdua de collita pot ser total. Però el dany no és únicament quantitatiu: el fong trenca la pell del raïm i genera compostos volàtils desagradables que <em class='term'>inhabiliten el fruit</em> per a l'elaboració de vins i caves DO. El cost dels tractaments preventius repetits, el combustible i el monitoratge suposen un increment notable dels costos d'explotació.",
    highlight: "El nostre biosensor permet passar de 'tractar per si de cas' a 'actuar on i quan hi ha infecció', optimitzant recursos i marges.",
    suggestions: ["resistència fungicides sofre","sostenibilitat agricultura precisió","oïdi erysiphe cendrosa"],
    confidence: 0.94
  },
  {
    id: 7,
    keys: ["clima","climàtic","temperatura","calor","hivern","supervivència","miceli","cicle","escalfament","primavera","infecció","global"],
    title: "Canvi Climàtic i Oïdi",
    body: "L'increment de temperatures mitjanes a Catalunya altera el cicle biològic d'<em>Erysiphe necator</em>. Els hiverns més suaus redueixen la mortalitat del miceli i de les cleistotècies, resultant en <b>infeccions primàries més precoces i agressives</b> a la primavera. A més, les onades de calor extremes limiten l'aplicació de sofre (fitotoxicitat), creant un doble risc: el fong és actiu però els tractaments estàndard no es poden aplicar.",
    highlight: "Paradoxa climàtica: quan més calor fa, menys podem tractar amb el producte ecològic principal, i el fong aprofita exactament aquesta finestra.",
    suggestions: ["resistència fungicides sofre","com funciona el biosensor","impacte econòmic collita vi"],
    confidence: 0.93
  },
  {
    id: 8,
    keys: ["sostenibilitat","ecologia","medi","ambient","optimització","precisió","ods","emissions","agricultura","preventius","recursos"],
    title: "Optimització i Sostenibilitat Agrícola",
    body: "La detecció precoç permet un canvi de paradigma: <b>agricultura de precisió</b>. En lloc de tractaments preventius sistemàtics, s'actua únicament on el sensor confirma la presència del fong. Això redueix la quantitat de fitosanitaris aplicats, el combustible de la maquinària, i les emissions de CO₂ associades. El projecte contribueix als ODS de la ONU, especialment ODS 2 (Fam Zero), ODS 12 (Producció i Consum Responsables) i ODS 13 (Acció Climàtica).",
    highlight: "Menys tractaments → menys costos → menys impacte ambiental → vi més net. Tothom hi guanya.",
    suggestions: ["impacte econòmic collita vi","canvi climàtic temperatura","resistència fungicides sofre"],
    confidence: 0.92
  },
  {
    id: 9,
    keys: ["igem","competició","concurs","biologia","sintètica","equip","urv","universitat","tarragona","wiki","synthetic","biology","projecte"],
    title: "Projecte iGEM URV 2025",
    body: "L'equip <b>iGEM URV</b> (Universitat Rovira i Virgili, Tarragona) participa a la competició internacional de biologia sintètica <em class='term'>iGEM 2025</em>. El nostre projecte aplica les eines de la biologia sintètica — disseny de proteïnes, enginyeria genètica i electrònica biosensora — per resoldre un problema real i local: la gestió de l'oïdi a la viticultura catalana. La wiki recull tota la documentació: laboratori, model, human practices i contribucions a la comunitat científica.",
    highlight: "La viticultura és el sector agrari més emblemàtic de Catalunya. Protegir-la amb tecnologia pròpia és un acte de responsabilitat científica i cultural.",
    suggestions: ["com funciona el biosensor","compostos volàtils fong","proteïnes OBP receptor"],
    confidence: 0.9
  },
  {
    id: 10,
    keys: ["transducció","senyal","electric","elèctric","electrodes","amperometria","transductor","resposta","impedància","potenciometria","bluetooth","lora","wireless","alerta"],
    title: "Transducció del Senyal Bioquímic",
    body: "La detecció OBP–compost volàtil es transdueix en un <b>senyal elèctric</b> mitjançant un transductor electroquímic. El canvi conformacional de l'OBP en unir-se al lligand altera la impedància de la superfície del sensor. Aquesta variació d'impedància es mesura per <em class='term'>amperometria</em> o <em class='term'>potenciometria</em>. El senyal analògic es digitalitza i es transmet via mòdul wireless (Bluetooth/LoRa) al sistema de gestió de la finca.",
    highlight: "Del camp al mòbil en temps real: el viticultor rep una alerta geolocalitzada indicant quin sector de la vinya presenta risc d'infecció.",
    suggestions: ["com funciona el biosensor","proteïnes OBP receptor","compostos volàtils fong"],
    confidence: 0.96
  },
  {
    id: 11,
    keys: ["enginyeria","disseny","dbtl","build","test","learn","iteració","clonació","expressió","construir","provar","aprendre","metodologia","cicle"],
    title: "Cicles d'Enginyeria (DBTL)",
    body: "El projecte segueix la metodologia d'enginyeria sintètica <b>DBTL</b> (Disseny – Construir – Provar – Aprendre). En la fase de <em class='term'>Disseny</em> seleccionem i modelitzem les OBPs candidats. En la fase de <em class='term'>Construir</em> expressem la proteïna en sistemes heteròlegs (<em>E. coli</em> o llevat). En la fase de <em class='term'>Provar</em> validem la unió al lligand per fluorescència (FRET) i ITC. Finalment, <em class='term'>Aprenem</em> de cada iteració per optimitzar l'afinitat i la selectivitat.",
    highlight: "Cada cicle DBTL ens acosta més a un sensor funcional en condicions de camp reals.",
    suggestions: ["proteïnes OBP receptor","compostos volàtils fong","com funciona el biosensor"],
    confidence: 0.94
  }
];

// ============================================================
// SEARCH — token matching sobre camp pla, sense Fuse per a frases
// ============================================================

// Paraules de parada catalanes a ignorar en la cerca
const STOP = new Set(["com","que","els","les","del","dels","una","uns","unes","per","amb","sobre","quins","quines","quin","quina","molt","però","és","el","la","un","de","i","a","en","no","un","ens","ho","hi","seu","seva","tot","tots","tota","totes","han","hem","heu","era","eren","ser","ha","fas","fan","fer","pot","poden","puc","pots","vol","vull","volen","cal"]);

// Prepara un text pla per a cada entrada
KB.forEach(e => {
  const bodyText = e.body.replace(/<[^>]+>/g,' ').replace(/\s+/g,' ');
  e._flat = (e.keys.join(' ') + ' ' + e.title + ' ' + bodyText).toLowerCase()
    .normalize('NFD').replace(/[\u0300-\u036f]/g,''); // sense accents per comparació
});

function normalize(s) {
  return s.toLowerCase().normalize('NFD').replace(/[\u0300-\u036f]/g,'');
}

function tokenize(text) {
  return normalize(text)
    .replace(/[.,?!;:'"«»()]/g,' ')
    .split(/\s+/)
    .filter(w => w.length >= 3 && !STOP.has(w));
}

function search(query) {
  const tokens = tokenize(query);
  if (!tokens.length) return [];
  return KB
    .map(e => {
      const hits = tokens.filter(t => e._flat.includes(t)).length;
      return { e, hits };
    })
    .filter(r => r.hits > 0)
    .sort((a,b) => b.hits - a.hits || b.e.confidence - a.e.confidence);
}

// ============================================================
// CHAT
// ============================================================
const messagesEl = document.getElementById('messages');
const inputEl    = document.getElementById('chatInput');

setTimeout(() => {
  appendBotMsg({
    title: "Benvingut/da al sistema expert OïdiSense",
    body: "Sóc l'agent expert del projecte iGEM URV 2025. Pots preguntar-me sobre l'<b>oïdi</b> (<em class='term'>Erysiphe necator</em>), el nostre <b>biosensor</b> basat en proteïnes OBP, la <b>viticultura catalana</b>, compostos volàtils, impacte econòmic, canvi climàtic o qualsevol aspecte del projecte.",
    highlight: "Prova: «Com funciona el biosensor?» o «Quins compostos detecta el sensor?»",
    suggestions: ["com funciona el biosensor","proteïnes OBP receptor","impacte econòmic collita vi"],
    confidence: 1.0
  });
}, 300);

function sendMessage() {
  const q = inputEl.value.trim();
  if (!q) return;
  appendUserMsg(q);
  inputEl.value = '';
  showTyping();
  setTimeout(() => {
    removeTyping();
    const results = search(q);
    if (!results.length) {
      appendBotMsg({
        title: "No he trobat una resposta",
        body: "No tinc informació específica sobre <em>«" + escHtml(q) + "»</em>. Prova a reformular o utilitza els botons de temes.",
        highlight: null,
        suggestions: ["com funciona el biosensor","oïdi erysiphe cendrosa","compostos volàtils fong"],
        confidence: 0
      });
      return;
    }
    const best = results[0].e;
    appendBotMsg({ title: best.title, body: best.body, highlight: best.highlight, suggestions: best.suggestions, confidence: best.confidence });
  }, 500 + Math.random() * 300);
}

function appendUserMsg(text) {
  const d = document.createElement('div');
  d.className = 'msg user';
  d.innerHTML = `<div class="msg-avatar">👤</div><div class="msg-bubble">${escHtml(text)}</div>`;
  messagesEl.appendChild(d);
  scrollBottom();
}

function appendBotMsg({ title, body, highlight, suggestions, confidence }) {
  const d = document.createElement('div');
  d.className = 'msg bot';
  const pct = Math.round(confidence * 100);
  const confHtml = confidence > 0
    ? `<div class="confidence"><span>Rellevància</span><div class="conf-bar"><div class="conf-fill" style="width:${pct}%"></div></div><span>${pct}%</span></div>`
    : '';
  const hlHtml = highlight ? `<div class="highlight">${highlight}</div>` : '';
  const sugHtml = (suggestions||[]).length
    ? `<div class="suggestions-row">${suggestions.map(s=>`<span class="suggestion-chip" onclick="askQuick('${s.replace(/'/g,"\\'")}')"> ${s}</span>`).join('')}</div>`
    : '';
  d.innerHTML = `<div class="msg-avatar">🍇</div><div class="msg-bubble"><strong>${title}</strong>${body}${hlHtml}${confHtml}${sugHtml}</div>`;
  messagesEl.appendChild(d);
  scrollBottom();
}

let typingEl = null;
function showTyping() {
  typingEl = document.createElement('div');
  typingEl.className = 'msg bot';
  typingEl.innerHTML = `<div class="msg-avatar">🍇</div><div class="typing-dots"><span></span><span></span><span></span></div>`;
  messagesEl.appendChild(typingEl);
  scrollBottom();
}
function removeTyping() { if (typingEl) { typingEl.remove(); typingEl = null; } }
function askQuick(q) { inputEl.value = q; sendMessage(); }
function scrollBottom() { messagesEl.scrollTop = messagesEl.scrollHeight; }
function escHtml(s) { return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
inputEl.addEventListener('keydown', e => { if (e.key==='Enter') { e.preventDefault(); sendMessage(); } });
