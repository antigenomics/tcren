// Shared controls for static molecular and table views.
const nav = document.getElementById('view-controls');
const theme = document.createElement('select');
for (const [value, label] of [['system','System'], ['light','Light'], ['dark','Dark']]) {
    theme.add(new Option(label, value));
}
const label = document.createElement('label');
label.append('Theme ', theme);
nav.append(label);
const preference = matchMedia('(prefers-color-scheme: dark)');
theme.value = localStorage.getItem('tcren-dashboard-theme') || 'system';
function applyTheme() {
    document.documentElement.dataset.theme = theme.value === 'system'
        ? (preference.matches ? 'dark' : 'light') : theme.value;
}
applyTheme();
preference.addEventListener('change', applyTheme);
theme.addEventListener('change', () => {
    localStorage.setItem('tcren-dashboard-theme', theme.value);
    applyTheme();
});

const search = document.getElementById("benchmark-search");
if (search) {
  const label = document.getElementById("benchmark-pairing");
  const epitope = document.getElementById("benchmark-epitope");
  const table = document.getElementById("benchmark-table");
  const rows = [...table.tBodies[0].rows];
  function filter() {
    const query = search.value.trim().toLowerCase();
    const visible = [0, 0];
    for (const row of rows) {
      row.hidden = !(row.textContent.toLowerCase().includes(query) &&
        (!label.value || row.dataset.pairing === label.value) &&
        (!epitope.value || row.dataset.epitope === epitope.value));
      if (!row.hidden) visible[Number(row.dataset.pairing)]++;
    }
    document.getElementById("benchmark-visible").textContent =
      `${visible[0] + visible[1]} of ${rows.length} rows · ${visible[1]} ${label.options[1].text} · ${visible[0]} ${label.options[2].text}`;
  }
  for (const header of table.querySelectorAll("[data-sort]")) {
    header.addEventListener("click", () => {
      const column = Number(header.dataset.sort);
      const current = header.closest("th");
      const ascending = current.getAttribute("aria-sort") === "descending";
      table.querySelectorAll("[aria-sort]").forEach(cell => cell.setAttribute("aria-sort", "none"));
      current.setAttribute("aria-sort", ascending ? "ascending" : "descending");
      rows.sort((first, second) => {
        const a = Number(first.cells[column].textContent), b = Number(second.cells[column].textContent);
        if (!Number.isFinite(a)) return Number.isFinite(b) ? 1 : 0;
        if (!Number.isFinite(b)) return -1;
        return ascending ? a - b : b - a;
      });
      table.tBodies[0].append(...rows);
    });
  }
  search.addEventListener("input", filter);
  label.addEventListener("change", filter);
  epitope.addEventListener("change", filter);
  filter();
}

const response = await fetch('navigation.json').catch(() => null);
if (response?.ok) {
    const {links = [], corpora = []} = await response.json();
    if (corpora.length) {
        const corpus = document.createElement('select');
        for (const item of corpora) {
            const option = new Option(item.label, item.href, item.current, item.current);
            corpus.add(option);
        }
        corpus.addEventListener('change', () => location.assign(corpus.value));
        const chooser = document.createElement('label');
        chooser.append('Corpus ', corpus);
        nav.prepend(chooser);
    }
    for (const item of links) {
        const link = document.createElement('a');
        link.textContent = item.label;
        if (item.current) link.setAttribute('aria-current', 'page');
        const sync = () => {
            const target = new URL(item.href, location.href);
            if (item.model) {
                const current = new URL(location.href);
                for (const key of ['group', 'id']) {
                    if (current.searchParams.has(key)) target.searchParams.set(key, current.searchParams.get(key));
                }
            }
            link.href = target.href;
        };
        sync();
        document.addEventListener('tcren:model', sync);
        nav.append(link);
    }
}
